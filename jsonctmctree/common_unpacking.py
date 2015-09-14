"""
Code that is common for unpacking json inputs.

This is shared between the code that computes log likelihood and derivatives
and the code that computes posterior expectations of things like
linear combinations of labeled transitions on edges.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
import networkx as nx

__all__ = [
        'SimpleError',
        'SimpleShapeError',
        'get_observables_info',
        'get_tree_info',
        'get_prior_info',
        'get_dwell_info',
        'get_root_request_info',
        ]


# Simple errors can go directly into the json error message
# without a stack trace.


class SimpleError(Exception):
    pass


class SimpleShapeError(SimpleError):
    pass


def get_observables_info(j_in, nnodes, state_space_shape):
    naxes = state_space_shape.shape[0]
    observable_nodes = np.array(j_in['observable_nodes'])
    observable_axes = np.array(j_in['observable_axes'])
    iid_observations = np.array(j_in['iid_observations'])

    # check ndims
    if len(observable_nodes.shape) != 1:
        raise SimpleShapeError(
                'expected the array of observable nodes '
                'to be one-dimensional')
    if len(observable_axes.shape) != 1:
        raise SimpleShapeError(
                'expected the array of observable axes '
                'to be one-dimensional')

    # check ndim of iid observations, handling the empty case
    nobservables = observable_axes.shape[0]
    if len(iid_observations.shape) != 2:
        raise SimpleShapeError(
                'expected the array of i.i.d. observations '
                'to be two-dimensional when observations are available')

    # check dtypes
    if nobservables:
        if not issubclass(observable_nodes.dtype.type, np.integer):
            raise SimpleError(
                    'expected observable_nodes to be an array of integers, '
                    'but found dtype %s' % observable_nodes.dtype)
        if not issubclass(observable_axes.dtype.type, np.integer):
            raise SimpleError(
                    'expected observable_axes to be an array of integers, '
                    'but found dtype %s' % observable_axes.dtype)

    # check conformant shapes
    if observable_nodes.shape[0] != observable_axes.shape[0]:
        raise SimpleShapeError(
                'The array of observable_nodes has length %d '
                'and the array of observable_axes has length %d '
                'but these are expected to be identical.' % (
                    observable_nodes.shape[0],
                    observable_axes.shape[0]))
    if nobservables:
        if observable_nodes.shape[0] != iid_observations.shape[1]:
            raise SimpleShapeError(
                    'The array of observable_nodes has length %d '
                    'and each iid observation has length %d '
                    'but these are expected to be identical.' % (
                        observable_nodes.shape[0],
                        iid_observations.shape[1]))

    # check contents
    if (
            np.any(observable_nodes < 0) or
            np.any(observable_axes < 0) or
            np.any(iid_observations < 0)):
        raise SimpleError(
                'The arrays of observable_nodes, observable_axes, '
                'and iid_observations must have non-negative integers')
    if nobservables:
        if np.any(nnodes <= observable_nodes):
            raise SimpleError(
                    'Each node index in the observable_nodes sequence should '
                    'be an integer between 0 and node_count-1 inclusive. '
                    'One or more of the observable node indices are too large')
        if np.any(naxes <= observable_axes):
            raise SimpleError(
                    'Each axis index in the observable_axes sequence '
                    'should be an integer between 0 and the number of axes '
                    'of the state space, minus one, inclusive. '
                    'One or more of the observable axis indices are too large')
        if np.any(state_space_shape[observable_axes] <= iid_observations):
            raise SimpleError(
                    'The entries of each iid observation '
                    'should be within the range of observed axis '
                    'of the state space ')

    return (
            observable_nodes,
            observable_axes,
            iid_observations)


def get_prior_info(j_in):
    return (
            np.array(j_in['prior_feasible_states']),
            np.array(j_in['prior_distribution'], dtype=float))


def get_root_request_info(j_in):
    """
    Optionally return per-site root state information.

    For each site, return a linear combination of posterior
    probabilities of states at the root.
    This could be used for EM, for example.

    """
    states = j_in.get('root_posterior_states', None)
    expect = j_in.get('root_posterior_expect', None)
    none_count = sum(1 for x in (states, expect) if x is None)
    if none_count not in (0, 2):
        raise SimpleError('expected neither or both of '
                'root_posterior_states and '
                'root_posterior_expect to be provided')
    if not none_count:
        return np.array(states), np.array(expect, dtype=float)
    else:
        return None, None


def get_dwell_info(j_in):
    """These inputs are optional.
    """
    states = j_in.get('dwell_states', None)
    expect = j_in.get('dwell_expect', None)
    none_count = sum(1 for x in (states, expect) if x is None)
    if none_count not in (0, 2):
        raise SimpleError('expected neither or both of dwell_states '
                'and dwell_expect to be provided')
    if not none_count:
        return np.array(states), np.array(expect, dtype=float)
    else:
        return None, None


def _check_tree_row_indices(row, node_count):
    """This is just for error checking.
    """
    nodes = set(range(node_count))
    unexpected_row_indices = list(set(row) - nodes)
    if unexpected_row_indices:
        raise SimpleError(
                'Found unexpected row indices in the tree definition. '
                'Because the provided node count is %d, '
                'the row indices are expected to be non-negative integers '
                'less than %d. '
                'But the following row indices were observed: %s' % (
                    node_count, node_count, unexpected_row_indices))


def _check_tree_col_indices(col, node_count):
    """This is just for error checking.
    """
    nodes = set(range(node_count))
    unexpected_col_indices = list(set(col) - nodes)
    if unexpected_col_indices:
        raise SimpleError(
                'Found unexpected col indices in the tree definition. '
                'Because the provided node count is %d, '
                'the col indices are expected to be non-negative integers '
                'less than %d. '
                'But the following col indices were observed: %s' % (
                    node_count, node_count, unexpected_col_indices))


def get_tree_info(j_in):
    node_count = j_in['node_count']
    process_count = j_in['process_count']
    tree = j_in['tree']
    nodes = set(range(node_count))
    row = tree['row']
    col = tree['col']
    rate = np.array(tree['rate'], dtype=float)
    process = np.array(tree['process'])
    _check_tree_row_indices(row, node_count)
    _check_tree_col_indices(col, node_count)
    negative_rates = rate[rate < 0]
    if negative_rates:
        raise SimpleError(
                'the edge-specific rate scaling factors '
                'should be non-negative')
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
