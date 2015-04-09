"""
Test some invariants due to reversibility.

"""
from __future__ import division, print_function, absolute_import

from itertools import permutations, product

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import scipy.linalg

from jsonctmctree import expect
from jsonctmctree import interface
from jsonctmctree.testutil import sample_time_reversible_rate_matrix


def _process_ex(Q, d, observable_node, debug=False):
    """
    Use the more advanced interface.

    """
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

    dwell_states = [[0, 0], [0, 1], [0, 2]]
    dwell_expect = [1, 1, 1]

    scene = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = state_space_shape,
            root_prior = dict(
                states = states,
                probabilities = d.tolist()),
            tree = dict(
                row_nodes = nodes[:-1],
                column_nodes = nodes[1:],
                edge_processes = [0]*nedges,
                edge_rate_scaling_factors = [0.2]*nedges,
                ),
            process_definitions = [dict(
                row_states = [i for i, j in state_pairs],
                column_states = [j for i, j in state_pairs],
                transition_rates = transition_rates,
                )],
            observed_data = dict(
                nodes = [observable_node],
                variables = [1],
                iid_observations = [
                    [0],
                    [2],
                    [1],
                    [0],
                    [1],
                    ]))

    dwell_request = dict(
            property = 'ddwdwel',
            state_reduction = dict(
                states = dwell_states,
                weights = dwell_expect))

    transition_request = dict(
            property = 'ddntran',
            transition_reduction = dict(
                row_states = [i for i, j in state_pairs],
                column_states = [j for i, j in state_pairs],
                weights = [1 for i, j in state_pairs]))

    j_in = dict(
        scene=scene,
        requests=[dwell_request, transition_request])

    return interface.process_json_in(j_in, debug=debug)


def _process(Q, d, observable_node, debug=False):
    """
    Use the obsolete interface.

    """
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

    return expect.process_json_in(j_in, debug=debug)


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
    e0_trans = j0['edge_expectations']
    e1_trans = j1['edge_expectations']
    for site in range(nsites):
        assert_allclose(e0_trans[site], list(reversed(e1_trans[site])))

    # For each site, the dwell proportion expectations should be be the same,
    # but in reverse order along the edges.
    e0_dwell = j0['edge_dwell']
    e1_dwell = j1['edge_dwell']
    for site in range(nsites):
        assert_allclose(e0_dwell[site], list(reversed(e1_dwell[site])))

    # Compare to the more advanced interface.

    j0_ex = _process_ex(Q, d, nodes[0])
    j1_ex = _process_ex(Q, d, nodes[-1])

    e0_ex_dwell, e0_ex_trans = j0_ex['responses']
    e1_ex_dwell, e1_ex_trans = j1_ex['responses']

    for site in range(nsites):
        assert_allclose(e0_ex_trans[site], list(reversed(e1_ex_trans[site])))
        assert_allclose(e0_ex_dwell[site], list(reversed(e1_ex_dwell[site])))

    for site in range(nsites):
        assert_allclose(e0_ex_trans[site], e0_trans[site])
        assert_allclose(e1_ex_trans[site], e1_trans[site])
        assert_allclose(e0_ex_dwell[site], e0_dwell[site])
        assert_allclose(e1_ex_dwell[site], e1_dwell[site])
