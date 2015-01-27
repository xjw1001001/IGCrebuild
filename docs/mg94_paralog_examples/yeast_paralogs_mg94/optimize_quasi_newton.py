"""
This is another attempt to write a generic interface for optimization.

Some of the parameters are assumed to correspond to
logarithms of branch lengths.
The other parameters are assumed to have been transformed to
an unconstrained vector in n-dimensional Euclidean space.

"""
from __future__ import print_function, division

import functools
import copy

import numpy as np
import scipy.optimize

__all__ = ['minimize']

# TODO automated EM steps for edge rate scaling factors.
# Inspect the process definitions to get the exit rates
# and the transition sparsity patterns.

def edge_rate_EM(scene, iterations):
    """
    Update edge rate scaling factors using EM.

    Note that edge rates that are zero will cause problems
    because there is no conditional information
    about these edge rate scaling factors.
    If the EM edge rate estimate is like the conditional expectation
    of the transition count on an edge divided by the
    expected transition opportunity,
    then this ratio will be 0/0 for those edges.

    Parameters
    ----------
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.

    Returns
    -------
    edge_rates : sequence of floats
        These will not have been log transformed.

    """
    # Copy the scene.
    scene = copy.deepcopy(scene)

    # For each process definition,
    # compute the total exit rate from each state.
    # and compute the transition sparsity pattern.
    state_to_exit_rate_dicts = []
    process_count = scene['process_count']
    for process in range(process_count):
        process_definition = scene['process_definitions'][process]
        row_states = process_definition['row_states']
        column_states = process_definition['column_states']
        transition_rates = process_definition['transition_rates']
        transition_count = len(transition_rates)
        state_to_exit_rate = defaultdict(float)
        for transition_index in range(transition_count):
            row_state = tuple(row_states[transition_index])
            state_to_exit_rate[row_state] += transition_rates[transition_index]
        state_to_exit_rate_dicts.append(state_to_exit_rate)

    # Each iteration, compute expectations.
    for i in range(iterations):
        pass

    # Return the edge rates from the most recent EM calculation.
    return edge_rates


def _mixed_gradient_objective(
        scene,
        get_process_definitions,
        get_root_prior,
        verbose,
        nP, nB, X):
    """

    Parameters
    ----------
    scene : jsonctmctree scene dict
        This dictionary aggregates the statistical model and the observed data.
    get_process_definitions : user-provided function f(P)
        Returns a list of process definition dicts given the global parameters.
        The global parameters use the unbounded transformation.
    get_root_prior : user-provided function f(P)
        Returns a root_prior dict given the global parameters.
    nP : integer
        Dimensionality of unbounded transformation of global parameters.
    nB : integer
        Number of edge-specific rate scaling factors.
    verbose : bool
        Extra information is printed if this is True.
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
    delta = 1e-7

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
    log_likelihood_request = dict(property = 'SNNLOGL')
    derivatives_request = dict(property = 'SDNLOGL')
    requests = [log_likelihood_request, derivatives_request]

    # Create the jsonctmctree input dict,
    # requesting the log likelihood and some derivatives.
    j_in = dict(
            scene = scene,
            requests = [log_likelihood_request, derivatives_request])

    # Compute the negative log likelihood
    # and the part of its gradient related to branch lengths.
    j_out = jsonctmctree.interface.process_json_in(j_in)
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
        j_out = jsonctmctree.interface.process_json_in(j_in)
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


def minimize(
        scene,
        get_process_definitions,
        get_root_prior,
        verbose,
        P0, B0):
    # Parameters are as in similarly named
    # parameters of _mixed_gradient_objective.
    P0 = np.asarray(P0)
    B0 = np.asarray(B0)
    nP = P0.shape[0]
    nB = B0.shape[0]
    X0 = np.concatenate((P0, B0))
    func_and_grad = functools.partial(
            scene, get_process_definitions, get_root_prior, verbose,
            nP, nB)
    result = scipy.optimize.minimize(func_and_grad, X0, method='L-BFGS-B')
    return result, result.x[:nP], result.x[-nB:]
