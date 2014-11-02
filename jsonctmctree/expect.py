"""
Linear combinations of labeled transition conditional expectations.

These conditional expectations are computed per-edge per-site.
Site weights are not accepted,
because per-site information is returned anyway.
Neither log likelihood nor its derivatives are returned.

For each process, an extra vector the same size as the rate vector is required.

"""
from __future__ import division, print_function, absolute_import

import sys

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
        get_subtree_likelihoods)


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


def get_node_to_marginal_distn(
        f,
        node_to_subtree_array, distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        debug=False):
    """

    """
    # Deduce counts.
    assert_equal(len(distn.shape), 1)
    nstates = distn.shape[0]
    nsites = iid_observations.shape[0]

    # Precompute some stuff.
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # The likelihoods downstream of nodes have already been precomputed.
    
    # Build the marginal distributions.
    node_to_marginal_distn = {}

    ordered_nodes = list(get_node_evaluation_order(T, root))
    for node in reversed(ordered_nodes):

        if debug:
            print('  node', node, '...', file=sys.stderr)

        if node == root:
            next_distn = node_to_subtree_array[root] * distn[:, np.newaxis]
            tail_node = root

            col_sums_recip = pseudo_reciprocal(next_distn.sum(axis=0))
            next_distn = next_distn * col_sums_recip
        else:
            # For non-root nodes the 'upstream edge' is of interest,
            # because we want to compute the weighted sum of expectations
            # of labeled transitions along the edge.
            edge = child_to_edge[node]
            head_node, tail_node = edge
            edge_process = edge_to_process[edge]
            edge_rate = edge_to_rate[edge]

            # Ingredients:
            # Marginal distribution at the head node of the edge.
            # Probability transition matrix along the edge.
            # Partial likelihood of the subtree.

            # Extract the marginal state distribution
            # at the 'head' of the edge.
            head_marginal_distn = node_to_marginal_distn[head_node]
            subtree_array = node_to_subtree_array[tail_node]

            # This vectorized implementation was worked out in
            # one of the test files in this module.
            A = head_marginal_distn * pseudo_reciprocal(
                    f[edge_process].expm_mul(edge_rate, subtree_array))
            B = f[edge_process].expm_rmul(edge_rate, A.T)
            next_distn = B.T * subtree_array

        node_to_marginal_distn[tail_node] = next_distn

    return node_to_marginal_distn


def pseudo_reciprocal(A):
    """
    Elementwise function 0->0, x->1/x

    """
    A_filled = np.where(A, A, 1)
    return np.where(A, np.reciprocal(A_filled), 0)


def get_edge_to_site_expectations(
        nsites, nstates,
        f, expm_frechet_objects, node_to_marginal_distn,
        node_to_subtree_array, distn,
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        debug=False):
    """

    """
    # Precompute some stuff.
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # The likelihoods downstream of nodes have already been precomputed.
    # The marginal distributions at nodes have also been computed.

    # Compute the expectations on edges.
    # Skip the root because it has no associated edge.
    edge_to_site_expectations = {}
    ordered_nodes = list(get_node_evaluation_order(T, root))
    for node in reversed(ordered_nodes[:-1]):

        if debug:
            print('  node', node, '...', file=sys.stderr)

        # For non-root nodes the 'upstream edge' is of interest,
        # because we want to compute the weighted sum of expectations
        # of labeled transitions along the edge.
        edge = child_to_edge[node]
        edge_process = edge_to_process[edge]
        edge_rate = edge_to_rate[edge]
        head_node, tail_node = edge

        # Extract the marginal state distribution
        # at the 'head' of the edge.
        head_marginal_distn = node_to_marginal_distn[head_node]
        subtree_array = node_to_subtree_array[tail_node]

        # Explicitly compute the transition probability matrix
        # as the matrix exponential of the transition rate matrix.
        # Also compute a thing that I'm calling the K matrix,
        # whose i, j entry is the weighted sum of the expected
        # number of transitions on the edge, given that the state of the
        # node at the 'head' of the edge is i and that the
        # state of the node at the 'tail' of the edge is j.
        obj = expm_frechet_objects[edge_process]
        P, K = obj.get_expm_and_frechet(edge_rate)

        #FIXME I am using brute force; do something intelligent instead!

        # First get the J matrix for each site -- this is the joint
        # distribution over states at the endpoints of the edge.
        # It is brutally inefficient to compute this separately for
        # each site!
        # Next, compute the hadamard product of this J matrix
        # with the K matrix computed using the frechet derivative
        # of the matrix exponential.
        # Finally, take the sum of the all entries in this
        # hadamard product as the expectation associated with the site.
        # This is so bad!

        # Take entrywise reciprocal of P, taking care to not divide by zero.
        P_recip = pseudo_reciprocal(P)
        site_expectations = np.empty(nsites, dtype=float)

        for site in range(nsites):

            if debug:
                print('    site', site, '...', file=sys.stderr)

            #TODO please replace with a vectorized version
            J = np.zeros((nstates, nstates))
            for i in range(nstates):
                d = P[i] * subtree_array[:, site]
                total = d.sum()
                if total:
                    d /= total
                else:
                    d = np.zeros_like(d)
                p = head_marginal_distn[i, site]
                J[i] = p * d

            # Report site-specific K for debugging.
            # The i, j entry of this matrix should give the
            # requested weighted sum of transition counts on the edge.
            K_mod = K * P_recip

            site_expectations[site] = (J * K_mod).sum()

        edge_to_site_expectations[edge] = site_expectations

    return edge_to_site_expectations


def process_json_in(j_in, debug=False):

    if debug:
        print('unpacking json input...', file=sys.stderr)

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
    processes_row, processes_col, processes_rate, processes_expect = info

    # Deduce some counts.
    nstates = np.prod(state_space_shape)
    nsites = iid_observations.shape[0]

    # Interpret the prior distribution.
    feas = np.ravel_multi_index(prior_feasible_states.T, state_space_shape)
    distn = np.zeros(nstates, dtype=float)
    np.put(distn, feas, prior_distribution)

    if debug:
        print('initializing expm objects...', file=sys.stderr)

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
        obj = ExplicitExpmFrechet(state_space_shape, row, col, rate, expect)
        expm_frechet_objects.append(obj)

    # Always store the likelihood arrays.
    # The only purpose of this function is to compute the labeled
    # edge expecatations, and this requires precomputing
    # the likelihood arrays.
    store_all_likelihood_arrays = True

    # Precompute conditional likelihood arrays per node.
    if debug:
        print('computing subtree likelihoods...', file=sys.stderr)
    node_to_subtree_array = get_subtree_likelihoods(
            f, store_all_likelihood_arrays,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            state_space_shape,
            observable_nodes,
            observable_axes,
            iid_observations)

    # Check the shape of the array.
    # Avoid copying a lot of huge arrays.
    for node in range(nnodes):
        assert_equal(node_to_subtree_array[node].shape, (nstates, nsites))

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_subtree_array[root]
    likelihoods = distn.dot(arr)
    feasibilities = (likelihoods > 0)

    # If the problem is infeasible, then return without computing expectations.
    if not np.all(feasibilities):
        print(node_to_subtree_array)
        j_out = dict(
                status = 'success',
                feasibility = False,
                edge_expectations = None)
        return j_out
    
    # Compute the marginal distribution.
    if debug:
        print('computing marginal distributions...', file=sys.stderr)
    node_to_marginal_distn = get_node_to_marginal_distn(
            f,
            node_to_subtree_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            state_space_shape,
            observable_nodes,
            observable_axes,
            iid_observations,
            debug=debug)

    # Check the shape of the array.
    # Avoid copying a lot of huge arrays.
    for node in range(nnodes):
        assert_equal(node_to_marginal_distn[node].shape, (nstates, nsites))

    # Compute expectations.
    if debug:
        print('computing site expectations...', file=sys.stderr)
    edge_to_site_expectations = get_edge_to_site_expectations(
            nsites, nstates,
            f, expm_frechet_objects, node_to_marginal_distn,
            node_to_subtree_array, distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            state_space_shape,
            observable_nodes,
            observable_axes,
            iid_observations,
            debug=debug)

    # Map expectations back to edge indices.
    # Note that this is per site per edge.
    # Take the transpose of this, so that the outer index
    # loops over sites.
    expectations_out = []
    for edge in edges:
        site_expectations = edge_to_site_expectations[edge]
        expectations_out.append(site_expectations)
    expectations_out = zip(*expectations_out)


    # Create the output in a format that json will like.
    j_out = dict(
            status = 'success',
            feasibility = True,
            edge_expectations = expectations_out)

    return j_out
