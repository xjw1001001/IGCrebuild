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
weighted reduction (tran_reduction) that is defined on multivariate
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

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal

from .expm_helpers import ActionExpm
from . import expect
from . import ll

def process_json_in(j_in):
    """
    The part of the input that is the same across requests is as follows.
    I'm bundling all of this stuff together and calling it a 'scenario'.

    'scenario' : {
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

    # Interpret the prior distribution by converting it to a dense array.
    nstates = np.prod(toplevel.scenario.state_space_shape)
    feas = np.ravel_multi_index(
            toplevel.scenario.root_prior.states.T,
            toplevel.scenario.state_space_shape)
    distn = np.zeros(nstates, dtype=float)
    np.put(distn, feas, prior_distribution)

    # Unpack stuff related to the tree and its edges.
    info = get_tree_info(j_in)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info

    # For each process, precompute the objects that are capable
    # of computing expm_mul and rate_mul for log likelihoods
    # and for its derivative with respect to edge-specific rates.
    f = []
    for p in scenario.process_definitions:
        obj = ActionExpm(
                toplevel.scenario.state_space_shape,
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
            toplevel.scenario.state_space_shape,
            toplevel.scenario.observation_data.nodes,
            toplevel.scenario.observation_data.variables,
            toplevel.scenario.observation_data.iid_observations)

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_array[root]
    likelihoods = distn.dot(arr)
    
    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    requested_derivative_edge_indices = set(requested_derivatives)
    ei_to_derivatives = get_edge_derivatives(
            f, requested_derivative_edge_indices,
            node_to_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            state_space_shape,
            observable_nodes,
            observable_axes,
            iid_observations)

    # Apply the prior distribution and take logs of the likelihoods.
    log_likelihoods = np.log(likelihoods)

    # Adjust for infeasibility.
    feasibilities = np.isfinite(log_likelihoods)
    log_likelihoods = np.where(feasibilities, log_likelihoods, 0)

    # Reduce the log likelihoods according to the site weights.
    if np.all(feasibilities):
        feasibility = True
        log_likelihood = np.dot(site_weights, log_likelihoods)

        # Process the derivatives.
        # They are currently in the form of derivatives of edge rates
        # with respect to the likelihood, but we want to convert them to
        # derivatives of log likelihood with respect to the logs
        # of edge-specific rate scaling factors.
        # According to calculus we can do it as follows.
        # Also reduce the derivative array according to the site weights.
        edge_to_rate = dict(edge_rate_pairs)
        ei_to_d = {}
        for ei, derivatives in ei_to_derivatives.items():
            edge = edges[ei]
            edge_rate = edge_to_rate[edge]
            d = site_weights.dot(derivatives / likelihoods)
            ei_to_d[ei] = float(d)

        # Map the derivatives back to a list whose entries
        # match the requested order of the indices.
        derivatives_out = [ei_to_d[ei] for ei in requested_derivatives]

    else:
        feasibility = False
        log_likelihood = 0
        derivatives_out = [0] * len(requested_derivatives)


    # Create the output in a format that json will like.
    j_out = dict(
            status = 'success',
            feasibility = feasibility,
            log_likelihood = log_likelihood,
            edge_derivatives = derivatives_out)

    return j_out
