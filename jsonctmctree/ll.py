"""
Implement log likelihoods for complicated models.

This interface takes some care about memory usage,
while allowing more subtlety in the representation of observed data,
and while allowing more flexibility in the representation of
inhomogeneity of the process across branches.

"""
from __future__ import print_function, division

import argparse
import json
import traceback
import sys

import networkx as nx
import numpy as np
from numpy.testing import assert_equal
from scipy.linalg import expm, eig, inv
from scipy.sparse.linalg import expm_multiply
from scipy.sparse import coo_matrix

from node_ordering import get_node_evaluation_order
from expm_helpers import PadeExpm, EigenExpm, ActionExpm



def get_site_weights(j_in):
    return np.array(j_in['site_weights'])


def get_requested_derivatives(j_in):
    return np.array(j_in['requested_derivatives'])


def get_observables_info(j_in):
    return (
            np.array(j_in['observable_nodes']),
            np.array(j_in['observable_axes']),
            np.array(j_in['iid_observations']))


def get_prior_info(j_in):
    return (
            np.array(j_in['prior_feasible_states']),
            np.array(j_in['prior_distribution'], dtype=float))


def get_processes_info(j_in):
    processes_row = []
    processes_col = []
    processes_rate = []
    for j_process in j_in['processes']:
        processes_row.append(j_process['row'])
        processes_col.append(j_process['col'])
        processes_rate.append(j_process['rate'])
    return (
            np.array(processes_row),
            np.array(processes_col),
            np.array(processes_rate, dtype=float))


def get_tree_info(j_in):
    node_count = j_in['node_count']
    process_count = j_in['process_count']
    tree = j_in['tree']
    nodes = set(range(node_count))
    row = tree['row']
    col = tree['col']
    rate = np.array(tree['rate'], dtype=float)
    process = np.array(tree['process'])
    if not (set(row) <= nodes):
        raise Exception('unexpected node')
    if not (set(col) <= nodes):
        raise Exception('unexpected node')
    T = nx.DiGraph()
    T.add_nodes_from(range(node_count))
    edges = zip(row, col)
    T.add_edges_from(edges)
    if len(T.edges()) != len(edges):
        raise Exception('the tree has an unexpected number of edges')
    if len(edges) + 1 != len(T):
        raise Exception('expected the number of edges to be one more '
                'than the number of nodes')
    in_degree = T.in_degree()
    roots = [n for n in nodes if in_degree[n] == 0]
    if len(roots) != 1:
        raise Exception('expected exactly one root')
    for i in range(node_count):
        T.in_degree()
    root = roots[0]
    edges = zip(row, col)
    edge_rate_pairs = zip(edges, rate)
    edge_process_pairs = zip(edges, process)
    return T, root, edges, edge_rate_pairs, edge_process_pairs


def create_indicator_array(
        node,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations):
    """
    Create the initial array indicating observations.

    """
    nsites, nobservables = iid_observations.shape
    state_space_ndim = len(state_space_shape)
    state_space_axes = range(state_space_ndim)

    # Initialize the active array, initially in a high dimensional shape.
    # This array is large; for data with many iid sites,
    # such active arrays dominate the memory usage of the program.
    obs_shape = (nsites, ) + tuple(state_space_shape)
    obs = np.ones(obs_shape, dtype=float)

    # For each observable associated with the node under consideration,
    # apply the observation mask across all iid sites.
    local_observables = np.flatnonzero(observable_nodes == node)
    for idx in local_observables:
        states = iid_observations[:, idx]
        axis = observable_axes[idx]
        k = state_space_shape[axis]
        projection_shape = [k if i == axis else 1 for i in state_space_axes]
        mask_shape = (nsites, ) + tuple(projection_shape)
        ident = np.identity(k)
        obs *= np.take(ident, states, axis=0).reshape(mask_shape)

    # Reshape the observation array to 2d.
    # First collapse the dimensionality of the state space to 1d,
    # then transpose the array so that it has shape (nstates, nsites),
    # to prepare for P.dot(obs) where P has shape (nstates, nstates).
    obs = obs.reshape((nsites, np.prod(state_space_shape))).T

    # Return the observation indicator array.
    return obs


def get_conditional_likelihoods(
        # functions to compute expm_mul, per process (NEW)
        f,
        # boolean flag for storing all edges arrays (NEW)
        store_all,
        # original args
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        ):
    """
    Recursively compute conditional likelihoods at the root.

    Attempt to order things intelligently to avoid using
    more memory than is necessary.

    The data provided by the caller gives us a sparse matrix
    of shape (nsites, nnodes, nstates).

    """
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # For the few nodes that are active at a given point in the traversal,
    # we track a 2d array of shape (nsites, nstates).
    node_to_array = {}
    for node in get_node_evaluation_order(T, root):

        # When a node is activated, its associated array
        # is initialized to its observational likelihood array.
        arr = create_indicator_array(
                node,
                state_space_shape,
                observable_nodes,
                observable_axes,
                iid_observations)

        # When an internal node is activated,
        # this newly activated observational array is elementwise multiplied
        # by each of the active arrays of the child nodes.
        # The new elementwise product becomes the array
        # associated with the activated internal node.
        #
        # If we did not care about saving the per-node arrays,
        # then we could inactivate the child nodes and delete
        # their associated arrays, but because we want to re-use the
        # per-node arrays for edge length gradients, we keep them.
        for child in T.successors(node):
            arr *= node_to_array[child]
            if not store_all:
                del node_to_array[child]

        # When any node that is not the root is activated,
        # the matrix product P.dot(A) replaces A,
        # where A is the active array and P is the matrix exponential
        # associated with the parent edge.
        if node != root:
            edge = child_to_edge[node]
            edge_rate = edge_to_rate[edge]
            edge_process = edge_to_process[edge]
            arr = f[edge_process].expm_mul(edge_rate, arr)

        # Associate the array with the current node.
        node_to_array[node] = arr

    # If we had been deleting arrays as they become unnecessary for
    # the log likelihood calculation, then we would have only
    # a single active array remaining at this point, corresponding to the root.
    # But if we are saving the arrays for gradient calculations,
    # then we have more left.
    actual_keys = set(node_to_array)
    if store_all:
        desired_keys = set(T)
    else:
        desired_keys = {root}
    assert_equal(actual_keys, desired_keys)

    # Return the map from node to array.
    return node_to_array


def get_edge_derivative(
        # functions to compute expm_mul and rate_mul, per process (NEW),
        f,
        # single derivative edge (pair of head, tail nodes not index) (NEW)
        derivative_edge,
        # map from node to array returned by (NEW)
        node_to_array,
        # original args
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        ):
    """

    """
    # Some preprocessing.
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # Unpack the edge of interest.
    derivative_head_node, derivative_tail_node = derivative_edge

    # For each edge use a new edge-specific array.
    # Trace the evaluation back to the root.
    node_to_deriv_array = {}

    # Iterate over nodes on the path from the derivative edge tail node
    # up to the root of the tree.
    node = derivative_tail_node
    while True:

        # When a node is activated, its associated array
        # is initialized to its observational likelihood array.
        arr = create_indicator_array(
                node,
                state_space_shape,
                observable_nodes,
                observable_axes,
                iid_observations)

        # If we are not analyzing the root node then determine
        # the characteristics of the edge upstream of the node.
        if node != root:
            edge = child_to_edge[node]
            edge_rate = edge_to_rate[edge]
            edge_process = edge_to_process[edge]

        # At the tail node, compute the adjusted array
        # for the edge of interest.
        # At non-tail nodes, including the root,
        # compute array using the same method as in the derivative-free
        # evaluation.
        # Arrays that are specific to the edge of interest
        # may be deleted after they are used.
        if node == derivative_tail_node:
            arr = node_to_array[node]
            arr = f[edge_process].rate_mul(edge_rate, arr)
        else:
            for child in T.successors(node):
                if child in node_to_deriv_array:
                    arr *= node_to_deriv_array[child]
                    del node_to_deriv_array[child]
                else:
                    arr *= node_to_array[child]
            if node != root:
                arr = f[edge_process].expm_mul(edge_rate, arr)

        # Associate the array with the current node.
        node_to_deriv_array[node] = arr

        # If the current node is the root node then we are done.
        # Otherwise move to the parent node.
        if node == root:
            break
        else:
            head_node, tail_node = edge
            node = head_node

    # Among the arrays specific to the analysis of the derivative
    # of the edge of interest, only the array at the root
    # should still be active at this point.
    assert_equal(set(node_to_deriv_array), {root})
    return node_to_deriv_array[root]


def get_edge_derivatives(
        # functions to compute expm_mul and rate_mul, per process (NEW),
        f,
        # set of requested edge indices for edge rate derivatives (NEW)
        requested_derivative_edge_indices,
        # map from node to array returned by (NEW)
        node_to_array,
        # per-site likelihoods (1d array and not log likelihoods) (NEW)
        likelihoods,
        # distribution at the root (NEW)
        distn,
        # original args
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        ):
    """
    Recursively compute conditional likelihoods at the root.

    Attempt to order things intelligently to avoid using
    more memory than is necessary.

    The data provided by the caller gives us a sparse matrix
    of shape (nsites, nnodes, nstates).

    """
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # Compute the likelihood derivative for each requested edge length.
    edge_index_to_derivatives = dict()
    for edge_index in requested_derivative_edge_indices:

        # Get the edge that corresponds to the edge index of interest.
        derivative_edge = edges[edge_index]

        # Compute the array at the root node.
        arr = get_edge_derivative(
                f,
                derivative_edge,
                node_to_array,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                state_space_shape,
                observable_nodes,
                observable_axes,
                iid_observations)

        # Apply the prior distribution to the array.
        # Now we have the derivatives of the likelihoods with respect
        # to the edge-specific rate scaling factors.
        # We want the derivatives of the log likelihoods with respect
        # to the log of the edge-specific rate scaling factors.
        derivatives = distn.dot(arr)
        edge_index_to_derivatives[edge_index] = derivatives

    # Return the map from edge index to edge-specific derivatives.
    return edge_index_to_derivatives


def process_json_in(j_in):

    # Unpack some sizes and shapes.
    nnodes = j_in['node_count']
    nprocesses = j_in['process_count']
    state_space_shape = np.array(j_in['state_space_shape'])

    # Unpack site weights. (NEW)
    site_weights = get_site_weights(j_in)

    # Unpack the indices of the edges for which the derivatives
    # of the likelihood are requested. (NEW)
    requested_derivatives = get_requested_derivatives(j_in)

    # Unpack stuff related to observables.
    info = get_observables_info(j_in)
    observable_nodes, observable_axes, iid_observations = info

    # Unpack stuff related to the prior distribution.
    info = get_prior_info(j_in)
    prior_feasible_states, prior_distribution = info

    # Unpack stuff related to the edge-specific processes.
    info = get_processes_info(j_in)
    processes_row, processes_col, processes_rate = info

    # Unpack stuff related to the tree and its edges.
    info = get_tree_info(j_in)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info

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
    for edge_process in range(nprocesses):
        row = processes_row[edge_process]
        col = processes_col[edge_process]
        rate = processes_rate[edge_process]
        obj = expm_klass(state_space_shape, row, col, rate)
        f.append(obj)

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
    requested_derivative_edge_indices = set(requested_derivatives)
    ei_to_derivatives = get_edge_derivatives(
            f, requested_derivative_edge_indices,
            node_to_array, likelihoods, distn,
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


def main(args):
    try:
        s_in = sys.stdin.read()
        j_in = json.loads(s_in)
    except Exception as e:
        if args.debug:
            raise
        return dict(
                status = 'error',
                message = 'json parsing error: ' + traceback.format_exc())
    try:
        return process_json_in(j_in)
    except Exception as e:
        if args.debug:
            raise
        return dict(
                status = 'error',
                message = 'processing error: ' + traceback.format_exc())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true')
    j_out = main(parser.parse_args())
    print(json.dumps(j_out))
