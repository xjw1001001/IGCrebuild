"""
Brute force tests of weighted expectations on edges.

"""
from __future__ import division, print_function, absolute_import

from itertools import product

import numpy as np
from numpy.testing import assert_allclose, assert_array_less

from scipy.linalg import expm
from jsonctmctree.expect import process_json_in


def _get_tiny_model_marginal_distns(nnodes, r):
    """
    Models evolution of a binary state along a path.

    Return the marginal distributions computed by brute force.

    """
    nedges = nnodes - 1
    nodes = range(nnodes + 1)
    nstates = 2

    Q = np.array([
        [-r, r],
        [ 0, 0]], dtype=float)
    P = expm(Q)
    W = np.zeros((nnodes, nstates))
    for assign in product(range(nstates), repeat=nnodes):
        assign = list(assign)
        # compute a likelihood
        if assign[0] != 0:
            continue
        if assign[-1] != 1:
            continue
        lk = 1
        for i, j in zip(assign[:-1], assign[1:]):
            lk *= P[i, j]
        for node, state in enumerate(assign):
            W[node, state] += lk

    row_sums = W.sum(axis=1)
    return W / row_sums[:, np.newaxis]


def _get_tiny_model_expectations_summary(nnodes, r):
    """
    Models evolution of a binary state along a path.

    """
    nedges = nnodes - 1
    nodes = range(nnodes)
    nstates = 2

    j_in = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = [2],
            prior_feasible_states = [[0]],
            prior_distribution = [1.0],
            dwell_states = [[1]],
            dwell_expect = [1],
            tree = dict(
                row = nodes[:-1],
                col = nodes[1:],
                process = [0]*nedges,
                rate = [1]*nedges,
                ),
            processes = [dict(
                row = [[0]],
                col = [[1]],
                rate = [r],
                expect = [1],
                )],
            observable_nodes = nodes[-1:],
            observable_axes = [0],
            iid_observations = [[1]])
    return process_json_in(j_in)


def test_tiny_model():
    """
    Models evolution of a binary state along a path.

    """
    nnodes = 4
    nedges = nnodes - 1
    nodes = range(nnodes)
    nstates = 2
    nsites = 1
    site = 0

    dwells = []

    for r in (0.1, 0.2):
        j_out = _get_tiny_model_expectations_summary(nnodes, r)
        expectations = j_out['edge_expectations'][site]
        M = _get_tiny_model_marginal_distns(nnodes, r)
        # compute theoretical probabilities
        x = np.arange(nnodes)
        probs = (1 - np.exp(-r * x)) / (1 - np.exp(-r * nedges))

        # Check that expected transition counts on edges of the path
        # correspond exactly to marginal distribution changes.
        for i, j in zip(nodes[:-1], nodes[1:]):
            theoretical_decrease = probs[j] - probs[i]
            desired_decrease = M[j][1] - M[i][1]
            actual_decrease = expectations[i]
            assert_allclose(actual_decrease, desired_decrease)
            assert_allclose(actual_decrease, theoretical_decrease)

        # For each edge, check the dwell proportion for the absorbing state.
        dwell = np.array(j_out['edge_dwell'][site])

        # Assert that the dwell times are between 0 and 1 and are increasing.
        assert_array_less(0, dwell)
        assert_array_less(dwell, 1)
        assert_array_less(0, dwell[1:] - dwell[:-1])

        dwells.append(dwell)

    # Check that the dwell time in the absorbing state is greater
    # on each edge for the process whose rate is larger.
    assert_array_less(dwell[0], dwell[1])


def _get_simple_model_expectations_summary(nnodes, r0, r1, e0, e1):
    nedges = nnodes - 1
    nodes = range(nnodes)

    j_in = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = [3],
            prior_feasible_states = [[0]],
            prior_distribution = [1.0],
            tree = dict(
                row = nodes[:-1],
                col = nodes[1:],
                process = [0]*nedges,
                rate = [1]*nedges,
                ),
            processes = [dict(
                row = [[0], [0]],
                col = [[1], [2]],
                rate = [r0, r1],
                expect = [e0, e1])],
            observable_nodes = [nnodes-1],
            observable_axes = [0],
            iid_observations = [
                [0],
                [1],
                [0],
                [2]])
    return process_json_in(j_in)

def test_simple_model():

    # Define a tree that is a path with four nodes and three edges.
    # Only the end state is observable.
    # 0 --(0)-- 1 --(1)-- 2 --(3)-- 3
    # The rate matrix will have only three states,
    # consisting of an undecided state and two decided states.
    # The rate towards the decided state 2 is twice as fast
    # as the rate towards the decided state 1.

    nnodes = 4
    nstates = 3
    nedges = nnodes - 1
    r0 = 0.1
    r1 = 0.2
    r = r0 + r1

    # Compute some theoretical probabilities.
    x = np.arange(nnodes)
    theoretical_probs = (1 - np.exp(-r * x)) / (1 - np.exp(-r * nedges))
    theoretical_diffs = theoretical_probs[1:] - theoretical_probs[:-1]

    # Check the combined model.
    # Validate expectations at each of the four sites.
    j_out = _get_simple_model_expectations_summary(nnodes, r0, r1, 1, 1)
    edge_expectations = j_out['edge_expectations']
    assert_allclose(edge_expectations[0], 0)
    assert_allclose(edge_expectations[1], theoretical_diffs)
    assert_allclose(edge_expectations[2], 0)
    assert_allclose(edge_expectations[3], theoretical_diffs)

    # Check weighted expectations for only the first transition type.
    # Validate expectations at each of the four sites.
    j_out = _get_simple_model_expectations_summary(nnodes, r0, r1, 1, 0)
    edge_expectations = j_out['edge_expectations']
    assert_allclose(edge_expectations[0], 0)
    assert_allclose(edge_expectations[1], theoretical_diffs)
    assert_allclose(edge_expectations[2], 0)
    assert_allclose(edge_expectations[3], 0)

    # Check weighted expectations for only the second transition type.
    # Validate expectations at each of the four sites.
    j_out = _get_simple_model_expectations_summary(nnodes, r0, r1, 0, 1)
    edge_expectations = j_out['edge_expectations']
    assert_allclose(edge_expectations[0], 0)
    assert_allclose(edge_expectations[1], 0)
    assert_allclose(edge_expectations[2], 0)
    assert_allclose(edge_expectations[3], theoretical_diffs)

    # Check brute force calculation of marginal distributions.
    Q = np.array([
        [-r, r0, r1],
        [0, 0, 0],
        [0, 0, 0]], dtype=float)
    P = expm(Q)
    W = np.zeros((nnodes, nstates))
    for assign in product(range(nstates), repeat=nnodes):
        assign = list(assign)
        # compute a likelihood
        if assign[0] != 0:
            continue
        if assign[-1] != 1:
            continue
        lk = 1
        for i, j in zip(assign[:-1], assign[1:]):
            lk *= P[i, j]
        for node, state in enumerate(assign):
            W[node, state] += lk
    distn = W / W.sum(axis=1)[:, np.newaxis]
    brute_diffs = (distn[1:] - distn[:-1])[:, 1]
    assert_allclose(brute_diffs, theoretical_diffs)
