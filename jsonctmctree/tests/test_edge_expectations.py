"""
Brute force test of weighted expectations on edges.

"""
from __future__ import division, print_function, absolute_import

from jsonctmctree.expect import process_json_in


def test_simple_model():

    # Define a tree that is a path with three nodes and two edges.
    # Only the end state is observable.
    # The second edge has twice the rate of the first edge.
    # 0 --(0)-- 1 --(1)-- 2
    # The rate matrix will have only three states,
    # consisting of an undecided state and two decided states.
    # The rate towards the decided state 2 is twice as fast
    # as the rate towards the decided state 1.

    j_in = dict(
            node_count = 3,
            process_count = 1,
            state_space_shape = [3],
            prior_feasible_states = [[0]],
            prior_distribution = [1.0],
            tree = dict(
                row = [0, 1],
                col = [1, 2],
                process = [0, 0],
                rate = [1.0, 2.0]),
            processes = [dict(
                row = [[0], [0]],
                col = [[1], [2]],
                rate = [1.0, 2.0],
                expect = [1.0, 1.0])],
            observable_nodes = [2],
            observable_axes = [0],
            iid_observations = [
                [0],
                [1],
                [2]])

    j_out = process_json_in(j_in)
    print(j_out)
