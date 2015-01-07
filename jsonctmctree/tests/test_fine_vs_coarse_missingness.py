"""
Test fine-grained vs. coarse-grained missing data.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_allclose, assert_equal,
        assert_array_less, assert_)

from jsonctmctree.interface import process_json_in


def get_partial_scene():
    # The observed data entry is missing from the returned scene.
    scene = dict(
            node_count = 3,
            process_count = 1,
            state_space_shape = [2, 2],
            tree = dict(
                row_nodes = [0, 0],
                column_nodes = [1, 2],
                edge_rate_scaling_factors = [0.1, 0.2],
                edge_processes = [0, 0]),
            root_prior = dict(
                states = [[0, 0], [0, 1], [1, 1], [1, 0]],
                probabilities = [0.1, 0.2, 0.3, 0.4]),
            process_definitions = [dict(
                row_states = [[0, 0], [0, 1], [1, 1], [1, 0]],
                column_states = [[0, 1], [1, 1], [1, 0], [0, 0]],
                transition_rates = [1, 2, 3, 4])])
    return scene


def test_coarse_vs_fine_grained_missingness():
    # Check that log likelihoods are the same in both representations.

    # Get the analysis results for the coarse-grained missingness.
    scene = get_partial_scene()
    coarse_observations = dict(
        nodes = [0, 0, 1],
        variables = [0, 1, 0],
        iid_observations = [
            [0, 0, 1],
            [0, 1, 0],
            [0, 0, 0],
            [1, 1, 1],
            [0, 0, 0],
            [1, 0, 1]])
    scene['observed_data'] = coarse_observations
    request = {'property' : 'DNNLOGL'}
    j_in = dict(scene = scene, requests = [request])
    j_out_coarse = process_json_in(j_in)

    # Get the analysis results for the fine-grained missingness.
    scene = get_partial_scene()
    fine_observations = dict(
        nodes = [0, 0, 1, 1],
        variables = [0, 1, 0, 1],
        iid_observations = [
            [0, 0, 1, -1],
            [0, 1, 0, -1],
            [0, 0, 0, -1],
            [1, 1, 1, -1],
            [0, 0, 0, -1],
            [1, 0, 1, -1]])
    scene['observed_data'] = fine_observations
    request = {'property' : 'DNNLOGL'}
    j_in = dict(scene = scene, requests = [request])
    j_out_fine = process_json_in(j_in)

    # Assert that the log likelihoods are the same.
    # Note that this uses exact equality comparison for floating point,
    # so if this test fails in the future then consider allowing
    # some epsilon of closeness.
    assert_equal(j_out_coarse, j_out_fine)
