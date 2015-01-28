"""
This is an attempt to write a generic interface for optimization.

Some of the parameters are assumed to correspond to
logarithms of edge rate scaling factors.
The other parameters are assumed to have been transformed to
an unconstrained vector in n-dimensional Euclidean space.

"""
from __future__ import division, print_function, absolute_import

import collections
import functools
import copy
import sys

import numpy as np
from numpy.testing import assert_equal, assert_
import scipy.optimize

from . import interface

__all__ = ['optimize_quasi_newton', 'optimize_em']



def _get_transition_reduction(process_definition):
    """
    Define a transition reduction.

    The reduction will simply count the expected transitions,
    without distinguishing between transition types.

    """
    transition_reduction = dict(
            row_states = process_definition['row_states'],
            column_states = process_definition['column_states'],
            weights = [1 for r in process_definition['transition_rates']])
    return transition_reduction


def _get_state_reduction(process_definition):
    """
    Define a state reduction.

    This will use the pre-scaling exit rates as weights.

    """
    row_states = process_definition['row_states']
    transition_rates = process_definition['transition_rates']
    transition_count = len(transition_rates)
    state_to_exit_rate = collections.defaultdict(float)
    for transition_index in range(transition_count):
        row_state = tuple(row_states[transition_index])
        state_to_exit_rate[row_state] += transition_rates[transition_index]
    states = []
    exit_rates = []
    for state, exit_rate in state_to_exit_rate.items():
        states.append(list(state))
        exit_rates.append(exit_rate)
    state_reduction = dict(
            states = states,
            weights = exit_rates)
    return state_reduction


def _do_em_iteration(
        scene,
        observation_reduction,
        transition_reductions,
        state_reductions):
    """
    For each unique edge_process, request summaries for all edges.

    Parameters
    ----------
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.
    observation_reduction : dict defining site-specific weights, or None.
        A reduction over observations, or None for an unweighted summation.
    transition_reductions : sequence of transition reduction dicts
        A transition reduction is provided for each unique edge process.
    state_reductions : sequence of state reduction dicts
        A state reduction is provided for each unique edge process.

    Returns
    -------
    edge_rates : sequence of floats
        The edge rate scaling factors calculated in this EM iteration.

    """
    process_count = scene['process_count']
    node_count = scene['node_count']
    edge_count = node_count - 1
    edge_rates = [None] * edge_count
    edge_processes = scene['tree']['edge_processes']
    for i in range(process_count):

        # Extract the transition reduction and the state reduction.
        # Use these to define expectation requests.
        if observation_reduction is not None:
            trans_request = dict(
                    property = 'WDNTRAN',
                    observation_reduction = observation_reduction,
                    transition_reduction = transition_reductions[i])
            dwell_request = dict(
                    property = 'WDWDWEL',
                    observation_reduction = observation_reduction,
                    state_reduction = state_reductions[i])
        else:
            trans_request = dict(
                    property = 'SDNTRAN',
                    transition_reduction = transition_reductions[i])
            dwell_request = dict(
                    property = 'SDWDWEL',
                    state_reduction = state_reductions[i])

        # Compute the expectations.
        j_in = dict(
                scene = scene,
                requests = [trans_request, dwell_request])
        j_out = interface.process_json_in(j_in)
        trans_response, dwell_response = j_out['responses']

        # Update the edge rates for edges whose associated process index is i.
        for edge_idx in range(edge_count):
            if edge_processes[edge_idx] == i:
                transitions = trans_response[edge_idx]
                opportunity = dwell_response[edge_idx]
                edge_rate = transitions / opportunity
                edge_rates[edge_idx] = edge_rate

    # Assert that every edge rate has been defined.
    for edge_rate in edge_rates:
        assert_(edge_rate is not None)

    # Return the new edge rates for this EM iteration.
    return edge_rates


def optimize_em(scene, observation_reduction, iterations):
    """
    Update edge rate scaling factors using EM.

    Note that edge rates that are zero will cause problems
    because there is no conditional information
    about these edge rate scaling factors.
    If the EM edge rate estimate is like the conditional expectation
    of the transition count on an edge divided by the
    expected transition opportunity,
    then this ratio will be 0/0 for those edges.

    Separate reductions are used per unique process.
    In other words if the process_count is k, then each iteration
    of this EM will compute k transition expectations
    and k dwell expectations.

    Parameters
    ----------
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.
    observation_reduction : dict defining site-specific weights, or None.
        A reduction over observations, or None for an unweighted summation.
    iterations : integer
        Do this many EM iterations.

    Returns
    -------
    edge_rates : sequence of floats
        These will not have been log transformed.

    """
    # Copy the scene.
    scene = copy.deepcopy(scene)
    process_count = scene['process_count']

    # For each unique process definition,
    # get the transition reduction and the state reduction.
    transition_reductions = []
    state_reductions = []
    for process_defn in scene['process_definitions']:
        transition_reductions.append(_get_transition_reduction(process_defn))
        state_reductions.append(_get_state_reduction(process_defn))

    # Update edge rates using a few EM iterations.
    edge_rates = None
    for em_iteration_idx in range(iterations):

        # If new edge rates are available then update the scene.
        if edge_rates is not None:
            scene['tree']['edge_rate_scaling_factors'] = edge_rates

        # Compute new edge rates with a single EM iteration.
        edge_rates = _do_em_iteration(
                scene,
                observation_reduction,
                transition_reductions,
                state_reductions,
                )

    # Return the edge rates from the most recent EM calculation.
    return edge_rates


def _mixed_gradient_objective(
        verbose,
        scene,
        observation_reduction,
        get_process_definitions,
        get_root_prior,
        nP, nB, X):
    """

    Parameters
    ----------
    verbose : bool
        Extra information is printed if this is True.
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.
    observation_reduction : dict defining site-specific weights, or None
        A reduction over observations, or None.
        If this is None then an unweighted summation will be used instead.
    get_process_definitions : user-provided function f(P)
        Returns a list of process definition dicts given the global parameters.
        The global parameters use the unbounded transformation.
    get_root_prior : user-provided function f(P)
        Returns a root_prior dict given the global parameters.
    nP : integer
        Dimensionality of unbounded transformation of global parameters.
    nB : integer
        Number of edge-specific rate scaling factors.
    X : 1-d array of floats
        The full parameter vector used by the quasi-Newton search.

    Returns
    -------
    y : float
        Value of the objective function.
        For example this could be the negative log likelihood.
    dydX : 1-d array of floats
        This is the gradient.
        It consists of derivatives of the objective function
        with respect to the transformed parameters and edge-specific
        rate scaling factors.

    """
    # Use this difference for finite differences.
    delta = 1e-8

    # Unpack the quasi-Newton search vector.
    X = np.asarray(X)
    assert_equal(X.shape, (nP + nB, ))
    P = X[:nP]
    B = X[-nB:]
    edge_rates = np.exp(B)

    # Copy the scene.
    # Update the process definitions.
    # Update the root distribution.
    # Update the edge rate scaling factors.
    scene = copy.deepcopy(scene)
    scene['process_definitions'] = get_process_definitions(P)
    scene['root_prior'] = get_root_prior(P)
    scene['tree']['edge_rate_scaling_factors'] = edge_rates

    # Define the log likelihood request and the gradient request.
    if observation_reduction is not None:
        log_likelihood_request = dict(
                property = 'WNNLOGL',
                observation_reduction = observation_reduction)
        derivatives_request = dict(
                property = 'WDNDERI',
                observation_reduction = observation_reduction)
    else:
        log_likelihood_request = dict(property = 'SNNLOGL')
        derivatives_request = dict(property = 'SDNDERI')
    requests = [log_likelihood_request, derivatives_request]

    # Create the jsonctmctree input dict,
    # requesting the log likelihood and some derivatives.
    j_in = dict(
            scene = scene,
            requests = [log_likelihood_request, derivatives_request])

    # Compute the negative log likelihood
    # and the part of its gradient related to branch lengths.
    j_out = interface.process_json_in(j_in)
    responses = j_out['responses']
    neg_log_likelihood = -responses[0]
    dydB = [-x for x in responses[1]]

    # For each non-edge-specific parameter,
    # numerically estimate a derivative using finite differences.
    dydP = []
    if verbose:
        print(
                'computing finite differences for tree-wide parameters...',
                file=sys.stderr)
    for i in range(nP):
        P2 = P.copy()
        P2[i] += delta
        scene['process_definitions'] = get_process_definitions(P2)
        scene['root_prior'] = get_root_prior(P2)
        j_in = dict(
                scene = scene,
                requests = [log_likelihood_request])
        j_out = interface.process_json_in(j_in)
        negative_log_likelihood_i = -j_out['responses'][0]
        deriv = (negative_log_likelihood_i - neg_log_likelihood) / delta
        dydP.append(deriv)
        if verbose:
            print(negative_log_likelihood_i, deriv, file=sys.stderr)

    # Concatenate the finite-differences derivates w.r.t. global parameters
    # and the more explicitly computed derivatives w.r.t. logs of
    # edge rate scaling factors.
    dydX = np.concatenate((dydP, dydB))

    # Return the objective function and its gradient.
    y = neg_log_likelihood
    return y, dydX


def optimize_quasi_newton(
        verbose,
        scene,
        observation_reduction,
        get_process_definitions,
        get_root_prior,
        P0, B0):
    """
    Use a quasi-Newton search.

    Parameters
    ----------
    verbose : bool
        Extra information is printed if this is True.
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.
    observation_reduction : dict defining site-specific weights, or None
        A reduction over observations, or None.
        If this is None then an unweighted summation will be used instead.
    get_process_definitions : user-provided function f(P)
        Returns a list of process definition dicts given the global parameters.
        The global parameters use the unbounded transformation.
    get_root_prior : user-provided function f(P)
        Returns a root_prior dict given the global parameters.
    nP : integer
        Dimensionality of unbounded transformation of global parameters.
    nB : integer
        Number of edge-specific rate scaling factors.
    X : 1-d array of floats
        The full parameter vector used by the quasi-Newton search.

    Returns
    -------
    result : OptimizeResult
        The output of the underlying scipy.optimize search.
    P : sequence of floats
        The transformed non-edge parameters.
    B : sequence of floats
        The logs of edge rate scaling factors.

    """
    P0 = np.asarray(P0)
    B0 = np.asarray(B0)
    nP = P0.shape[0]
    nB = B0.shape[0]
    X0 = np.concatenate((P0, B0))
    func_and_grad = functools.partial(
            _mixed_gradient_objective,
            verbose,
            scene,
            observation_reduction,
            get_process_definitions,
            get_root_prior,
            nP, nB)
    result = scipy.optimize.minimize(
            func_and_grad, X0, jac=True, method='L-BFGS-B')
    return result, result.x[:nP], result.x[-nB:]
