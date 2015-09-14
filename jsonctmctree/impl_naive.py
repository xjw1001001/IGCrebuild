"""
Naive implementation of interface.py.

This naive implementation eagerly precomputes the full array associated
with all base properties, and then subsequently returns reductions
requested by the caller.  The naiveity comes into play by not cleverly
incorporating the reductions into the calculations.
Although this implementation is not efficient, it should be useful for
testing the correctness of more efficient implementations.

"""
from __future__ import division, print_function, absolute_import

import sys

import numpy as np
from numpy.testing import assert_equal

from .expm_helpers import (
        ActionExpm,
        ImplicitDwellExpmFrechet,
        ImplicitTransitionExpmFrechetEx,
        )
from .common_likelihood import (
        get_conditional_likelihoods, get_subtree_likelihoods)
from .common_unpacking_ex import TopLevel, interpret_tree, interpret_root_prior
from .common_reduction import apply_reductions
from . import expect
from . import ll


def _eagerly_precompute_dwell_objects(scene):
    """
    Precompute dwell times for each state on each edge at each site.

    This will be done more cleverly in the less naive implementation later.
    Predefine the dwell objects for each process for each site.

    Returns a nested list so that arr[i][j] is the dwell object
    for integer state i and associated with process j.

    """
    nstates = np.prod(scene.state_space_shape)

    dwell_objects = []
    for dwell_state_index in range(nstates):

        # The linear combination of states is not very interesting.
        arr = []
        dwell_state = np.array(np.unravel_index(
                dwell_state_index, scene.state_space_shape))
        dwell_states = np.array([dwell_state])
        dwell_weights = np.ones(1)

        # For this dwell state index track one object per unique process.
        for p in scene.process_definitions:
            obj = ImplicitDwellExpmFrechet(
                    scene.state_space_shape,
                    p.row_states,
                    p.column_states,
                    p.transition_rates,
                    dwell_states,
                    dwell_weights,
                    )
            arr.append(obj)
        dwell_objects.append(arr)

    return dwell_objects


def _apply_eagerly_precomputed_dwell_objects(
        scene,
        expm_objects, dwell_objects, node_to_marginal_distn,
        node_to_subtree_likelihoods, prior_distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        ):
    """

    """
    nprocesses = len(scene.process_definitions)
    nsites = scene.observed_data.iid_observations.shape[0]
    nstates = np.prod(scene.state_space_shape)
    assert_equal(len(dwell_objects), nprocesses)

    edge_to_dwell_expectations = expect.get_edge_to_site_expectations(
            nsites, nstates,
            expm_objects, dwell_objects, node_to_marginal_distn,
            node_to_subtree_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            scene.state_space_shape,
            scene.observed_data.nodes,
            scene.observed_data.variables,
            scene.observed_data.iid_observations,
            debug=False)

    # These dwell times will be scaled by the edge-specific scaling factor.
    # We want to remove that effect.
    for edge, edge_rate in edge_rate_pairs:
        if edge_rate:
            edge_to_dwell_expectations[edge] /= edge_rate

    # Map expectations back to edge indices.
    # The output will be like (nedges, nsites).
    edge_dwell_out = []
    for edge in edges:
        dwell_expectations = edge_to_dwell_expectations[edge]
        edge_dwell_out.append(dwell_expectations)

    return edge_dwell_out


def _compute_transition_expectations(
        scene,
        expm_objects, expm_frechet_objects, node_to_marginal_distn,
        node_to_subtree_likelihoods, prior_distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        debug=False,
        ):
    """

    """
    nprocesses = len(scene.process_definitions)
    nsites = scene.observed_data.iid_observations.shape[0]
    nstates = np.prod(scene.state_space_shape)

    edge_to_site_expectations = expect.get_edge_to_site_expectations(
            nsites, nstates,
            expm_objects, expm_frechet_objects, node_to_marginal_distn,
            node_to_subtree_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            scene.state_space_shape,
            scene.observed_data.nodes,
            scene.observed_data.variables,
            scene.observed_data.iid_observations,
            debug=debug)

    # Map expectations back to edge indices.
    # This will have shape (nedges, nsites).
    expectations_out = []
    for edge in edges:
        site_expectations = edge_to_site_expectations[edge]
        expectations_out.append(site_expectations)

    return expectations_out


def process_json_in(j_in, debug=False):
    """
    Eagerly computes everything and then later applies the reductions.

    A less naive method would avoid precomputed unnecessary things
    and would possibly move some reductions inside the computations.

    """
    toplevel = TopLevel(j_in)

    # Precompute the size of the state space.
    nstates = np.prod(toplevel.scene.state_space_shape)
    iid_observation_count = (
            toplevel.scene.observed_data.iid_observations.shape[0])

    # Interpret the prior distribution by converting it to a dense array.
    prior_distn = interpret_root_prior(toplevel.scene)

    # Interpret stuff related to the tree and its edges.
    info = interpret_tree(toplevel.scene)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info
    nedges = len(edges)

    # For each process, precompute the objects that are capable
    # of computing expm_mul and rate_mul for log likelihoods
    # and for its derivative with respect to edge-specific rates.
    expm_objects = []
    for p in toplevel.scene.process_definitions:
        obj = ActionExpm(
                toplevel.scene.state_space_shape,
                p.row_states,
                p.column_states,
                p.transition_rates)
        expm_objects.append(obj)

    # In this naive implementation, always store all intermediate arrays.
    store_all = True
    node_to_subtree_likelihoods = get_subtree_likelihoods(
            expm_objects,
            store_all,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations,
            )

    # In this naive implementation, always store all intermediate arrays.
    store_all = True
    node_to_conditional_likelihoods = get_conditional_likelihoods(
            expm_objects,
            store_all,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations,
            )

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_conditional_likelihoods[root]
    likelihoods = prior_distn.dot(arr)
    assert_equal(len(likelihoods.shape), 1)

    # If the likelihood at any site is zero
    # then report infeasibility for the scene.
    if not np.all(likelihoods):
        return dict(
                status = 'infeasible',
                responses = None)

    # Apply the prior distribution and take logs of the likelihoods.
    log_likelihoods = np.log(likelihoods)

    # Compute the marginal distributions.
    # The full node array will have shape (nsites, nstatates, nnodes).
    if debug:
        print('computing marginal distributions...', file=sys.stderr)
    node_to_marginal_distn = expect.get_node_to_marginal_distn(
            expm_objects,
            node_to_subtree_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations,
            debug=debug)
    arr = []
    for i in range(toplevel.scene.node_count):
        arr.append(node_to_marginal_distn[i])
    full_node_array = np.array(arr).T

    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    # In this naive implementation compute all edge derivatives.
    requested_derivative_edge_indices = set(range(nedges))
    ei_to_derivatives = ll.get_edge_derivatives(
            expm_objects, requested_derivative_edge_indices,
            node_to_conditional_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Fill an array with all unreduced derivatives.
    derivatives = np.empty((iid_observation_count, nedges))
    for ei, der in ei_to_derivatives.items():
        derivatives[:, ei] = der / likelihoods

    # Precompute dwell objects without regard to the requests.
    # This will obviously be changed for the less naive implementation.
    all_dwell_objects = _eagerly_precompute_dwell_objects(toplevel.scene)

    # Apply the precomputed dwell objects,
    # creating an array like (nsites, nedges, nstates) after the transposition.
    # A less naive implementation would compute this less eagerly.
    arr = []
    for dwell_state_index in range(nstates):
        dwell_objects = all_dwell_objects[dwell_state_index]
        arr.append(_apply_eagerly_precomputed_dwell_objects(
                toplevel.scene,
                expm_objects, dwell_objects, node_to_marginal_distn,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs))
    full_dwell_array = np.array(arr).T

    # Respond to each request, using a 'scene' common to all requests.
    responses = []
    for req in toplevel.requests:
        prefix, suffix = req.property[:3], req.property[-4:]
        observation_code, edge_code, state_code = prefix

        if suffix == 'logl':
            out = log_likelihoods
        elif suffix == 'deri':
            out = derivatives
        elif suffix == 'dwel':
            out = full_dwell_array
        elif suffix == 'tran':
            expm_transition_objects = []
            for p in toplevel.scene.process_definitions:
                obj = ImplicitTransitionExpmFrechetEx(
                        toplevel.scene.state_space_shape,
                        p.row_states,
                        p.column_states,
                        p.transition_rates,
                        req.transition_reduction.row_states,
                        req.transition_reduction.column_states,
                        req.transition_reduction.weights,
                        )
                expm_transition_objects.append(obj)
            arr = _compute_transition_expectations(
                toplevel.scene,
                expm_objects, expm_transition_objects, node_to_marginal_distn,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                debug=debug)
            out = np.array(arr).T
        elif suffix == 'root':
            out = full_node_array[:, :, root]
        elif suffix == 'node':
            out = full_node_array

        # Apply reductions along axes according to the request.
        out = apply_reductions(toplevel.scene.state_space_shape, req, out)

        # Convert the ndarray to a list.
        # If the response is a zero ndim ndarray, then this will
        # convert to a single number.
        # Otherwise if the response is not an ndarray
        # (for example if it is a Python float),
        # then an AttributeError will be raised and the response
        # will be unchanged.
        try:
            out = out.tolist()
        except AttributeError as e:
            pass

        # Append the response.
        responses.append(out)

    # Return the response.
    return dict(status='feasible', responses=responses)
