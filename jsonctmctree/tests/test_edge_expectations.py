"""
Brute force test of weighted expectations on edges.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

from jsonctmctree.expect import process_json_in


def test_tiny_model():

    # Define a tree that is a path.
    # The model transitions from an initial state to an absorbing state.
    nnodes = 4
    nedges = nnodes - 1
    
    for r in (0.1, 0.2):
        j_in = dict(
                node_count = nnodes,
                process_count = 1,
                state_space_shape = [2],
                prior_feasible_states = [[0]],
                prior_distribution = [1.0],
                tree = dict(
                    row = np.arange(nedges).tolist(),
                    col = np.arange(1, nedges+1).tolist(),
                    process = np.zeros(nedges, dtype=int).tolist(),
                    rate = np.ones(nedges, dtype=int).tolist(),
                    ),
                processes = [dict(
                    row = [[0]],
                    col = [[1]],
                    rate = [r],
                    expect = [1],
                    )],
                observable_nodes = [nnodes-1],
                observable_axes = [0],
                iid_observations = [[1]])

        j_out = process_json_in(j_in)
        print(j_out)


def test_simple_model():

    # Define a tree that is a path with four nodes and three edges.
    # Only the end state is observable.
    # 0 --(0)-- 1 --(1)-- 2 --(3)-- 3
    # The rate matrix will have only three states,
    # consisting of an undecided state and two decided states.
    # The rate towards the decided state 2 is twice as fast
    # as the rate towards the decided state 1.

    nnodes = 4
    nedges = nnodes - 1

    j_in = dict(
            node_count = 4,
            process_count = 1,
            state_space_shape = [3],
            prior_feasible_states = [[0]],
            prior_distribution = [1.0],
            tree = dict(
                row = np.arange(nedges).tolist(),
                col = np.arange(1, nedges+1).tolist(),
                process = np.zeros(nedges, dtype=int).tolist(),
                rate = np.ones(nedges, dtype=int).tolist(),
                ),
            processes = [dict(
                row = [[0], [0]],
                col = [[1], [2]],
                rate = [0.1, 0.2],
                expect = [1.0, 1.0])],
            observable_nodes = [nnodes-1],
            observable_axes = [0],
            iid_observations = [
                [0],
                [1],
                [0],
                [2]])

    j_out = process_json_in(j_in)
    print(j_out)
