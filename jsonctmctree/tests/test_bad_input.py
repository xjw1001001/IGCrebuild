"""
Test input errors.

"""
from __future__ import print_function, division

from numpy.testing import assert_raises

import jsonctmctree.ll
from jsonctmctree.ll import SimpleError, SimpleShapeError


def get_good_univariate_input():
    return dict(
            node_count = 2,
            process_count = 1,
            state_space_shape = [2],
            prior_feasible_states = [[0], [1]],
            prior_distribution = [0.4, 0.6],
            tree = dict(
                process = [0],
                col = [0],
                row = [1],
                rate = [4.2]),
            requested_derivatives = [],
            processes = [dict(
                row = [[0], [1]],
                col = [[1], [0]],
                rate = [0.1, 10.0])],
            observable_nodes = [0, 1],
            observable_axes = [0, 0],
            site_weights = [1, 1, 1],
            iid_observations = [
                [0, 0],
                [0, 0],
                [1, 0]])


def test_good_univariate_input():
    j_in = get_good_univariate_input()
    j_out = jsonctmctree.ll.process_json_in(j_in)


def test_bad_observation_shape():
    j_in = get_good_univariate_input()
    j_in['iid_observations'] = [0, 0, 1]
    assert_raises(SimpleShapeError, jsonctmctree.ll.process_json_in, j_in)


def test_out_of_bounds_observable_axis():
    j_in = get_good_univariate_input()
    j_in['observable_axes'] = [0, 1]
    assert_raises(SimpleError, jsonctmctree.ll.process_json_in, j_in)
