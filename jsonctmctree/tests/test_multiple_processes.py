"""
Test some properties of bivariate processes on different edges.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_allclose, assert_equal,
        assert_array_less, assert_)

from jsonctmctree.ll import process_json_in


def _get_json_input(process_indices):
    assert_(process_indices in ((0, 0), (0, 1), (1, 0), (1, 1)))
    state_space_shape = (2, 2)
    nstates = np.prod(state_space_shape)
    nnodes = 3
    nedges = nnodes - 1
    nodes = range(nnodes)
    edges = range(nedges)
    prior_feasible_states = [[0, 0], [1, 1]]
    prior_distribution = [0.5, 0.5]

    # Define two distinct bivariate processes.
    # The first bivariate process has independent components,
    # so this first process is really just like two univariate process.
    # The second bivariate process has correlated components.
    fast = 2.0
    slow = 1.0
    processes = [
            dict(
                row = [
                    [0, 0], [0, 0],
                    [0, 1], [0, 1],
                    [1, 0], [1, 0],
                    [1, 1], [1, 1],
                    ],
                col = [
                    [0, 1], [1, 0],
                    [0, 0], [1, 1],
                    [0, 0], [1, 1],
                    [0, 1], [1, 0],
                    ],
                rate = [
                    slow, slow,
                    slow, slow,
                    slow, slow,
                    slow, slow,
                    ],
                ),
            dict(
                row = [
                    [0, 0], [0, 0],
                    [0, 1], [0, 1],
                    [1, 0], [1, 0],
                    [1, 1], [1, 1],
                    ],
                col = [
                    [0, 1], [1, 0],
                    [0, 0], [1, 1],
                    [0, 0], [1, 1],
                    [0, 1], [1, 0],
                    ],
                rate = [
                    slow, slow,
                    fast, fast,
                    fast, fast,
                    slow, slow,
                    ],
                )
            ]

    # Nothing is observable at the root.
    # Only one variable of the bivariate process is observable
    # at the first leaf.
    # Both variables of the bivariate process are observable
    # at the second leaf.
    observable_nodes = [1, 2, 2]
    observable_axes = [0, 0, 1]
    iid_observations = [
            [0, 0, 0],
            [1, 1, 1],
            [1, 1, 1],
            [0, 0, 0],
            ]

    j_in = dict(
            node_count = nnodes,
            process_count = 2,
            state_space_shape = state_space_shape,
            prior_feasible_states = prior_feasible_states,
            prior_distribution = prior_distribution,
            tree = dict(
                row = [0, 0],
                col = [1, 2],
                process = process_indices,
                rate = [0.25]*nedges,
                ),
            processes = processes,
            observable_nodes = observable_nodes,
            observable_axes = observable_axes,
            iid_observations = iid_observations,
            site_weights = [1]*len(iid_observations),
            requested_derivatives = [0, 1],
            )
    return j_in


def test_multiple_processes():
    ll00 = process_json_in(_get_json_input((0, 0)))['log_likelihood']
    ll01 = process_json_in(_get_json_input((0, 1)))['log_likelihood']
    ll10 = process_json_in(_get_json_input((1, 0)))['log_likelihood']
    ll11 = process_json_in(_get_json_input((1, 1)))['log_likelihood']
    assert_allclose(ll00, ll10)
    assert_allclose(ll01, ll11)
    assert_array_less(ll00, ll11)
