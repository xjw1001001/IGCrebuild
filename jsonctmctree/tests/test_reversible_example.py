"""
Test some invariants due to reversibility.

"""
from __future__ import division, print_function, absolute_import

from itertools import permutations, product

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import scipy.linalg

from jsonctmctree.expect import process_json_in
from jsonctmctree.sampling import(
        assert_symmetric_matrix,
        sample_symmetric_rates,
        sample_distn,
        sample_time_reversible_rate_matrix,
        sample_time_nonreversible_rate_matrix,
        )


def _process(Q, d, observable_node):
    state_space_shape = (2, 3)
    nstates = np.prod(state_space_shape)
    nnodes = 4
    nedges = nnodes - 1
    nodes = range(nnodes)
    edges = range(nedges)
    states = list(product(
        range(state_space_shape[0]),
        range(state_space_shape[1]),
        ))
    state_pairs = list(permutations(states, 2))
    ntrans = len(state_pairs)
    assert_equal(ntrans, nstates * (nstates - 1))
    row, col = zip(*state_pairs)
    idx_row = np.ravel_multi_index(np.transpose(row), state_space_shape)
    idx_col = np.ravel_multi_index(np.transpose(col), state_space_shape)
    transition_rates = [Q[i, j] for i, j in zip(idx_row, idx_col)]

    root_posterior_states = [[0, 0]]
    root_posterior_expect = [1]

    dwell_states = [[0, 0], [0, 1], [0, 2]]
    dwell_expect = [1, 1, 1]

    j_in = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = state_space_shape,
            prior_feasible_states = states,
            prior_distribution = d.tolist(),
            dwell_states = dwell_states,
            dwell_expect = dwell_expect,
            root_posterior_states = root_posterior_states,
            root_posterior_expect = root_posterior_expect,
            tree = dict(
                row = nodes[:-1],
                col = nodes[1:],
                process = [0]*nedges,
                rate = [0.2]*nedges,
                ),
            processes = [dict(
                row = [i for i, j in state_pairs],
                col = [j for i, j in state_pairs],
                rate = transition_rates,
                expect = [1]*ntrans,
                )],
            observable_nodes = [observable_node],
            observable_axes = [1],
            iid_observations = [
                [0],
                [2],
                [1],
                [0],
                [1],
                ])
    return process_json_in(j_in)


def test_time_reversible_invariants():
    # Use an interesting state space shape but a less interesting tree (a path).
    #
    # Use a time-reversible matrix, and compare posterior statistics
    # when information is known at only the final node in the chain vs.
    # when information is known at only the initial node in the chain.
    #
    nsites = 5
    state_space_shape = np.array([2, 3])
    nstates = np.prod(state_space_shape)
    nnodes = 4
    nedges = nnodes - 1
    nodes = range(nnodes)
    edges = range(nedges)
    Q, d = sample_time_reversible_rate_matrix(nstates)

    # First look at posterior statistics when the final node is known
    # but the initial node is set to the equilibrium distribution.
    # Then look at posterior statistics when then final node state is known.
    j0 = _process(Q, d, nodes[0])
    j1 = _process(Q, d, nodes[-1])

    # For each site, the transition count expectations should be be the same,
    # but in reverse order along the edges.
    e0 = j0['edge_expectations']
    e1 = j1['edge_expectations']
    for site in range(nsites):
        assert_allclose(e0[site], list(reversed(e1[site])))

    # For each site, the dwell proportion expectations should be be the same,
    # but in reverse order along the edges.
    e0 = j0['edge_dwell']
    e1 = j1['edge_dwell']
    for site in range(nsites):
        assert_allclose(e0[site], list(reversed(e1[site])))

