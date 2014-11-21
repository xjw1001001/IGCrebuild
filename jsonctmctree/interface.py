"""
Begin a new interface.

This includes support for user requests for six 'posterior' base properties:
    * LOGL: log likelihood
    * DERI: derivatives with respect to log edge rates
    * TRAN: transition count expectations
    * DWEL: dwell proportion expectations
    * ROOT: state count expectations at the root
    * NODE: state count expectations at all nodes

Each base property is extended to allow one or more reductions:
    * reduction across iid observations (observation_reduction)
    * reduction across edges (edge_reduction)
    * reduction across states (state_reduction)

Each requested reduction may be defined by either simple summation
or by a user-specified weighted sum over a user-specified sequence of indices.

The transition base property request always requires an additional
weighted reduction (transition_reduction) that is defined on multivariate
state pairs.

For each property request,
the response array will have one axis for each D in the property prefix.
If the array has zero axes then a single floating point number will be returned.

The extended properties have 7-letter names according to the following scheme:
observation (1) | edge (1) | state (1) | base property name (4)
where the letters of the 3-letter prefix are from:
(D)istinct (no reduction)
(S)um (unweighted sum)
(W)eighted sum
(N)ot applicable

The 6 base properties can be extended as follows to a total of 39 properties:
{D,S,W}NNLOGL : 3
{D,S,W}DNDERI : 3
{D,S,W}{D,W}{D,W}DWEL : 12
{D,S,W}{D,S,W}NTRAN : 9
{D,S,W}N{D,W}ROOT : 6
{D,S,W}N{D,W}NODE : 6

The interface is limited in that it does not support the following:
    * continuous observations along time intervals
    * non-axis-aligned state aggregate observations
    * noisy observations (subsumes state aggregate observations)
    * partitions of observations within a single scene

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal

from .expm_helpers import (
        ActionExpm,
        ImplicitDwellExpmFrechet,
        ImplicitTransitionExpmFrechetEx,
        )
from .common_likelihood import get_conditional_likelihoods
from .common_unpacking_ex import TopLevel, interpret_tree, interpret_root_prior
from . import expect
from . import ll


def _sparse_reduction(A, indices, weights, axis):
    """
    A weighted sum along an axis.

    Parameters
    ----------
    A : ndarray
        the ndarray to reduce
    indices : 1d ndarray
        the indices to be used in the reduction along the axis of interest
    weights : 1d ndarray
        the array of weights to be used in the reduction
    axis : int
        the axis along which to reduce the input ndarray

    """
    A = np.asarray(A)
    indices = np.asarray(indices)
    weights = np.asarray(weights)
    assert_equal(indices.ndim, 1)
    assert_equal(weights.ndim, 1)
    assert_equal(indices.shape, weights.shape)
    return np.tensordot(
            np.take(A, indices, axis=axis),
            weights,
            axes=([axis], [0]))


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
        node_to_subtree_array, distn,
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
            node_to_subtree_array, distn,
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
        node_to_subtree_array, distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        ):
    """

    """
    nprocesses = len(scene.process_definitions)
    nsites = scene.observed_data.iid_observations.shape[0]
    nstates = np.prod(scene.state_space_shape)

    edge_to_site_expectations = expect.get_edge_to_site_expectations(
            nsites, nstates,
            expm_objects, expm_frechet_objects, node_to_marginal_distn,
            node_to_subtree_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            scene.state_space_shape,
            scene.observed_data.nodes,
            scene.observed_data.variables,
            scene.observed_data.iid_observations,
            debug=False)

    # Map expectations back to edge indices.
    # This will have shape (nedges, nsites).
    expectations_out = []
    for edge in edges:
        site_expectations = edge_to_site_expectations[edge]
        expectations_out.append(site_expectations)

    return expectations_out


def _process_json_in_naive(j_in, debug=False):
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
    distn = interpret_root_prior(toplevel.scene)

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

    # Precompute conditional likelihood arrays per node.
    # In this naive implementation, always store all intermediate arrays.
    store_all_likelihood_arrays = True
    node_to_array = get_conditional_likelihoods(
            expm_objects, store_all_likelihood_arrays,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_array[root]
    likelihoods = distn.dot(arr)
    assert_equal(len(likelihoods.shape), 1)

    # If the likelihood at any site is zero
    # then report infeasibility for the scene.
    if not np.all(likelihoods):
        return dict(
                status = 'infeasible',
                responses = None)

    # Apply the prior distribution and take logs of the likelihoods.
    log_likelihoods = np.log(likelihoods)

    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    # In this naive implementation compute all edge derivatives.
    requested_derivative_edge_indices = set(range(nedges))
    ei_to_derivatives = ll.get_edge_derivatives(
            expm_objects, requested_derivative_edge_indices,
            node_to_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Fill an array with all unreduced derivatives.
    derivatives = np.empty((iid_observation_count, nedges))
    for ei, der in ei_to_derivatives.items():
        derivatives[:, ei] = der / likelihoods

    # Compute the marginal distributions.
    # The full node array will have shape (nsites, nstatates, nnodes).
    if debug:
        print('computing marginal distributions...', file=sys.stderr)
    node_to_marginal_distn = expect.get_node_to_marginal_distn(
            expm_objects,
            node_to_array, distn,
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

    # Precompute dwell objects without regard to the requests.
    # This will obviously be changed for the less naive implementation.
    all_dwell_objects = _eagerly_precompute_dwell_objects(toplevel.scene)

    # Eagerly apply the eagerly computed dwell objects,
    # creating an array like (nsites, nedges, nstates) after the transposition.
    # A less naive implementation would compute this less eagerly.
    arr = []
    for dwell_state_index in range(nstates):
        dwell_objects = all_dwell_objects[dwell_state_index]
        arr.append(_apply_eagerly_precomputed_dwell_objects(
                toplevel.scene,
                expm_objects, dwell_objects, node_to_marginal_distn,
                node_to_array, distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs))
    full_dwell_array = np.array(arr).T

    #{D,S,W}NNLOGL : 3
    #{D,S,W}DNDERI : 3
    #{D,S,W}{D,W}{D,W}DWEL : 12
    #{D,S,W}{D,S,W}NTRAN : 9
    #{D,S,W}N{D,W}ROOT : 6
    #{D,S,W}N{D,W}NODE : 6

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
                node_to_array, distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs)
            out = np.array(arr).T
        elif suffix == 'root':
            out = full_node_array[:, :, root]
        elif suffix == 'node':
            out = full_node_array

        # Prepare to apply the reductions.
        reduction_axis = 0

        # Apply the observation reduction if any.
        if observation_code == 'd':
            reduction_axis += 1
        elif observation_code == 's':
            out = np.sum(out, axis=reduction_axis)
        elif observation_code == 'w':
            indices = req.observation_reduction.observation_indices
            weights = req.observation_reduction.weights
            out = _sparse_reduction(out, indices, weights, reduction_axis)

        # Apply the edge reduction if any.
        if edge_code == 'd':
            reduction_axis += 1
        elif edge_code == 's':
            out = np.sum(out, axis=reduction_axis)
        elif edge_code == 'w':
            indices = req.edge_reduction.edges
            weights = req.edge_reduction.weights
            out = _sparse_reduction(out, indices, weights, reduction_axis)

        # Apply the state reduction if any.
        if state_code == 'd':
            reduction_axis += 1
        elif state_code == 's':
            out = np.sum(out, axis=reduction_axis)
        elif state_code == 'w':
            indices = np.ravel_multi_index(
                    req.state_reduction.states.T,
                    toplevel.scene.state_space_shape)
            weights = req.state_reduction.weights
            out = _sparse_reduction(out, indices, weights, reduction_axis)

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


def process_json_in(j_in):
    """
    The part of the input that is the same across requests is as follows.
    I'm bundling all of this stuff together and calling it a 'scene'.

    'scene' : {
        'node_count' : 4,
        'process_count' : 2,
        'state_space_shape' : [2, 2],
        'root_prior' : {
            'states' : ...,
            'probabilities' : ...},
        'tree' : {
            'row_nodes' : ...,
            'column_nodes' : ...,
            'edge_rate_scaling_factors' : ...,
            'edge_processes' : ...},
        'process_definitions' : [
            {
                'row_states' : ...,
                'column_states' : ...,
                'transition_rates' : ...}, ...],
        'observed_data' : {
            'nodes' : [...],
            'variables' : [...],
            'iid_observations' : [[...]]}
    }

    The requests part of the input is an array of json objects,
    each of which has a 'property' (one of the 39 properties listed above)
    and may have one or more weighted reduction members.
    The number of weighted reduction definition members is equal to the
    number of 'w' characters in the 3-letter prefix of the property.
    Each weighted reduction definition is an object with an array
    of indices and an array of weights.

    'requests' : [
        {
            'property' : 'snnlogl'
        },
        {
            'property' : 'wwwdwel',
            'observation_reduction' : {
                'observation_indices' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'edge_reduction' : {
                'edges' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'state_reduction' : {
                'states' : [[0, 0], [0, 1], [0, 2]],
                'weights' : [1, 1, 1]}
        },
        {
            'property' : 'wsntran',
            'observation_reduction' : {
                'observation_indices' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'transition_reduction' : {
                'row_states' : [[0, 0], [0, 1], [0, 2]],
                'column_states' : [[0, 1], [0, 2], [0, 0]],
                'weights' : [1, 1, 1]}
        }
        ]

    """
    return _process_json_in_naive(j_in, debug=False)

    # Convert the parsed json input into something that looks
    # more like an object, with some input validation and some type converion.
    toplevel = TopLevel(j_in)
    """
        'root_prior' : {
            'states' : ...,
            'probabilities' : ...},
        'tree' : {
            'row_nodes' : ...,
            'column_nodes' : ...,
            'edge_rate_scaling_factors' : ...,
            'edge_processes' : ...},
        'process_definitions' : [
            {
                'row_states' : ...,
                'column_states' : ...,
                'transition_rates' : ...}, ...],
        'observed_data' : {

    """
    # Precompute the size of the state space.
    nstates = np.prod(toplevel.scene.state_space_shape)
    iid_observation_count = (
            toplevel.scene.observed_data.iid_observations.shape[0])

    # Interpret the prior distribution by converting it to a dense array.
    distn = interpret_root_prior(toplevel.scene)

    # Interpret stuff related to the tree and its edges.
    info = interpret_tree(toplevel.scene)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info
    nedges = len(edges)

    # For each process, precompute the objects that are capable
    # of computing expm_mul and rate_mul for log likelihoods
    # and for its derivative with respect to edge-specific rates.
    f = []
    for p in toplevel.scene.process_definitions:
        obj = ActionExpm(
                toplevel.scene.state_space_shape,
                p.row_states,
                p.column_states,
                p.transition_rates)
        f.append(obj)

    # Determine whether to store intermediate arrays
    # or whether to store only as many as necessary to compute
    # the log likelihood.
    #store_all_likelihood_arrays = bool(np.any(requested_derivatives))
    #FIXME
    store_all_likelihood_arrays = True

    # Precompute conditional likelihood arrays per node.
    node_to_array = get_conditional_likelihoods(
            f, store_all_likelihood_arrays,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_array[root]
    likelihoods = distn.dot(arr)
    assert_equal(len(likelihoods.shape), 1)

    # If the likelihood at any site is zero
    # then report infeasibility for the scene.
    if not np.all(likelihoods):
        return dict(
                status = 'infeasible',
                responses = None)

    # Determine which edge derivatives are requested.
    """
    requested_derivative_edge_indices = set()
    for r in toplevel.requests:
        if r.property.endswith('deri'):
            requested_derivative_edge_indices.update(r.edges)
    """
    # FIXME for now we just requre all edge indices
    requested_derivative_edge_indices = set(range(nedges))
    
    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    ei_to_derivatives = ll.get_edge_derivatives(
            f, requested_derivative_edge_indices,
            node_to_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Apply the prior distribution and take logs of the likelihoods.
    log_likelihoods = np.log(likelihoods)

    # Reduce the log likelihoods according to the site weights.
    #log_likelihood = np.dot(site_weights, log_likelihoods)

    # Process the derivatives.
    # They are currently in the form of derivatives of edge rates
    # with respect to the likelihood, but we want to convert them to
    # derivatives of log likelihood with respect to the logs
    # of edge-specific rate scaling factors.
    # According to calculus we can do it as follows.
    # Also reduce the derivative array according to the site weights.
    #ei_to_d = {}
    #for ei, derivatives in ei_to_derivatives.items():
        #edge = edges[ei]
        #d = site_weights.dot(derivatives / likelihoods)
        #ei_to_d[ei] = float(d)

    # Fill an array with unreduced derivatives.
    derivatives = np.empty((iid_observation_count, nedges))
    for ei, der in ei_to_derivatives.items():
        derivatives[:, ei] = der / likelihoods

    # Map the derivatives back to a list whose entries
    # match the requested order of the indices.
    #derivatives_out = [ei_to_d[ei] for ei in requested_derivatives]

    #{D,S,W}NNLOGL : 3
    #{D,S,W}DNDERI : 3
    #{D,S,W}{D,W}{D,W}DWEL : 12
    #{D,S,W}{D,S,W}NTRAN : 9
    #{D,S,W}N{D,W}ROOT : 6
    #{D,S,W}N{D,W}NODE : 6

    # Create one response for each request.
    responses = []
    for req in toplevel.requests:
        prefix, suffix = req.property[:3], req.property[-4:]
        observation_code, edge_code, state_code = prefix
        if suffix == 'logl':
            response_per_site = log_likelihoods
        elif suffix == 'deri':
            response_per_site = derivatives
        elif suffix == 'dwel':
            #FIXME
            response_per_site = log_likelihoods
        elif suffix == 'tran':
            #FIXME
            response_per_site = log_likelihoods
        elif suffix == 'root':
            #FIXME
            response_per_site = log_likelihoods
        elif suffix == 'node':
            #FIXME
            response_per_site = log_likelihoods

        if observation_code == 'd':
            response = response_per_site
        elif observation_code == 's':
            response = np.sum(response_per_site, axis=0)
        elif observation_code == 'w':
            indices = req.observation_reduction.observation_indices
            weights = req.observation_reduction.weights
            arr = np.take(response_per_site, indices, axis=0)
            response = np.dot(arr, weights)

        responses.append(response)

    # Return the response.
    return dict(status='feasible', responses=responses)
