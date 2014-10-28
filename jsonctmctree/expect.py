"""
Linear combinations of labeled transition conditional expectations.

These conditional expectations are computed per-edge per-site.
Site weights are not accepted,
because per-site information is returned anyway.
Neither log likelihood nor its derivatives are returned.

For each process, an extra vector the same size as the rate vector is required.

"""
from __future__ import division, print_function, absolute_import

import networkx as nx
import numpy as np
from numpy.testing import assert_equal

from .expm_helpers import (
        PadeExpm, EigenExpm, ActionExpm,
        ExplicitExpmFrechet)

from .expm_helpers import create_dense_rate_matrix

from .node_ordering import get_node_evaluation_order

from .common_unpacking import (
        SimpleError,
        SimpleShapeError,
        get_observables_info,
        get_tree_info,
        get_prior_info)

from .common_likelihood import (
        create_indicator_array,
        get_conditional_likelihoods)


#TODO remove
class ExplicitExpmFrechet(object):
    """
    This is for computing conditional expectations on edges.

    """
    def __init__(self, state_space_shape, row, col, rate, coeffs):
        self.Q = create_sparse_rate_matrix(state_space_shape, row, col, rate)
        self.E = create_sparse_rate_matrix(state_space_shape, row, col, coeffs)

    def get_expm_and_frechet(self, rate_scaling_factor):
        QA = self.Q.A
        EA = self.E.A
        t = rate_scaling_factor
        P, interact = scipy.linalg.expm_frechet(QA * t, EA)
        del QA
        del EA
        return P, interact


def get_processes_info(j_in):
    processes_row = []
    processes_col = []
    processes_rate = []
    processes_expect = []
    for j_process in j_in['processes']:
        processes_row.append(j_process['row'])
        processes_col.append(j_process['col'])
        processes_rate.append(j_process['rate'])
        processes_expect.append(j_process['expect'])
    return (
            np.array(processes_row),
            np.array(processes_col),
            np.array(processes_rate, dtype=float),
            np.array(processes_expect, dtype=float))


def get_expectations(
        f, expm_frechet_objects,
        node_to_array, distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations):
    """

    """
    # Precompute some stuff.
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # The conditional likelihoods have already been precomputed.

    node_to_marginal_distn = {}

    ordered_nodes = list(get_node_evaluation_order(T, root))
    for node in reversed(ordered_nodes):
        upstream_edge_array = node_to_array[node]
        if node == root:
            likelihoods = distn.dot(upstream_edge_array)
            # It is possible that some likelihoods that should be zero
            # are listed as tiny floating point numbers,
            # so force these to be zero.
            # Then normalize the likelihoods as the marginal distribution
            # at the root.
            d = np.clip(likelihoods, 0, np.inf)
            node_to_marginal_distn[root] = d / d.sum()
        else:
            # For non-root nodes the 'upstream edge' is of interest,
            # because we want to compute the weighted sum of expectations
            # of labeled transitions along the edge.
            edge = child_to_edge[node]
            edge_process = edge_to_process[edge]
            edge_rate = edge_to_rate[edge]

            # Explicitly compute the transition probability matrix
            # as the matrix exponential of the transition rate matrix,
            # and also compute the 'interaction matrix' as its
            # frechet derivative in the direction of the provided
            # transition expectation coefficients.
            obj = expm_frechet_objects[edge_process]
            P, numerator = obj.get_expm_and_frechet(edge_rate)


    # Compute Frechet derivatives of matrix exponentials on-demand,
    # and be careful with memory so delete them when no longer needed.



def process_json_in(j_in):

    # Unpack some sizes and shapes.
    nnodes = j_in['node_count']
    nprocesses = j_in['process_count']
    state_space_shape = np.array(j_in['state_space_shape'])

    # Unpack stuff related to the tree and its edges.
    info = get_tree_info(j_in)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info

    # Unpack stuff related to observables.
    info = get_observables_info(j_in, nnodes, state_space_shape)
    observable_nodes, observable_axes, iid_observations = info

    # Unpack stuff related to the prior distribution.
    info = get_prior_info(j_in)
    prior_feasible_states, prior_distribution = info

    # Unpack stuff related to the edge-specific processes.
    info = get_processes_info(j_in)
    processes_row, processes_col, processes_rate = info

    # Interpret the prior distribution.
    nstates = np.prod(state_space_shape)
    feas = np.ravel_multi_index(prior_feasible_states.T, state_space_shape)
    distn = np.zeros(nstates, dtype=float)
    np.put(distn, feas, prior_distribution)

    # For each process, precompute the objects that are capable
    # of computing expm_mul and rate_mul for log likelihoods
    # and for its derivative with respect to edge-specific rates.
    #expm_klass = EigenExpm # TODO soft-code this
    #expm_klass = PadeExpm # TODO soft-code this
    expm_klass = ActionExpm # TODO soft-code this
    f = []
    expm_frechet_objects = []
    for edge_process in range(nprocesses):
        row = processes_row[edge_process]
        col = processes_col[edge_process]
        rate = processes_rate[edge_process]
        expect = processes_expect[edge_process]
        obj = expm_klass(state_space_shape, row, col, rate)
        f.append(obj)
        obj = ExplicitExpmFrechet(state_space_shape, row, col, expect)
        expm_frechet_objects.append(obj)

    # Determine whether to store intermediate arrays
    # or whether to store only as many as necessary to compute
    # the log likelihood.
    store_all_likelihood_arrays = bool(np.any(requested_derivatives))

    # Precompute conditional likelihood arrays per node.
    node_to_array = get_conditional_likelihoods(
            f, store_all_likelihood_arrays,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            state_space_shape,
            observable_nodes,
            observable_axes,
            iid_observations)

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_array[root]
    likelihoods = distn.dot(arr)
    
    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    #requested_derivative_edge_indices = set(requested_derivatives)
    #ei_to_derivatives = get_edge_derivatives(
    foo = get_expectations(
            f, expm_frechet_objects,
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
