"""
"""
from __future__ import division, print_function, absolute_import

from itertools import permutations

import numpy as np
from numpy.testing import assert_allclose, dec

import jsonctmctree.interface


def get_poisson_definition(nstates, poisson_rate):
    row_states = []
    column_states = []
    transition_rates = []
    for i, j in permutations(range(nstates), 2):
        row_states.append([i])
        column_states.append([j])
        rate = poisson_rate / nstates
        transition_rates.append(rate)
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates,
            )
    return process_definition


def get_star_tree(nleaves):
    nedges = nleaves
    root_node = nleaves
    tree = dict(
            row_nodes = [root_node] * nedges,
            column_nodes = range(nleaves),
            edge_rate_scaling_factors = [1] * nedges,
            edge_processes = [0] * nedges,
            )
    return tree


def get_informative_root_prior():
    root_prior = dict(
            states = [[0]],
            probabilities = [1],
            )
    return root_prior


def get_uniform_root_prior(nstates):
    root_prior = dict(
            states = [[i] for i in range(nstates)],
            probabilities = [1/nstates] * nstates,
            )
    return root_prior


def get_observed_data(nleaves):
    observed_data = dict(
            nodes = range(nleaves),
            variables = [0] * nleaves,
            iid_observations = [[1] * nleaves],
            )
    return observed_data


def get_poisson_scene(nleaves, nstates, poisson_rate):
    # nleaves : the number of leaves in the star tree
    # nstates : the size of the state space
    # poisson_rate : the rate of the poisson process
    #
    tree = get_star_tree(nleaves)
    root_prior = get_informative_root_prior()
    process_definition = get_poisson_definition(nstates, poisson_rate)
    observed_data = get_observed_data(nleaves)
    nnodes = nleaves + 1
    scene = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = [nstates],
            tree = tree,
            root_prior = root_prior,
            process_definitions = [process_definition],
            observed_data = observed_data,
            )
    return scene


def _check_poisson(nleaves, nstates, poisson_rate):
    desired_ll_per_leaf = np.log(-np.expm1(-poisson_rate) / nstates)
    desired_ll = nleaves * desired_ll_per_leaf
    scene = get_poisson_scene(nleaves, nstates, poisson_rate)
    request = dict(property = 'snnlogl')
    j_in = dict(scene = scene, requests = [request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    actual_ll = j_out['responses'][0]
    assert_allclose(actual_ll, desired_ll)


def test_numerically_easy_poisson():
    # This example should not require paying attention to scaling.
    nstates = 4
    poisson_rate = 1e-4
    nleaves = 10
    _check_poisson(nleaves, nstates, poisson_rate)


@dec.knownfailureif(True)
def test_numerically_difficult_poisson():
    # This example may require the implementation to pay attention to scaling.
    nstates = 4
    poisson_rate = 1e-4
    nleaves = 100
    _check_poisson(nleaves, nstates, poisson_rate)
