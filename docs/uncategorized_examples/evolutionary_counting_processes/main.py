"""
This is really a K80 model because the nucleotide distribution is uniform.

Reproduce an example from Minin and Suchard
"Counting labeled transitions in continuous-time Markov models of evolution."

"""
from __future__ import print_function, division

import itertools

import numpy as np
from numpy.testing import assert_allclose

from jsonctmctree.interface import process_json_in


def gen_K80():
    # Use the nucleotide order from the Minin and Suchard paper.
    # A, G, C, T
    transitions = ((0, 1), (1, 0), (2, 3), (3, 2))
    for i in range(4):
        for j in range(4):
            if i != j:
                if (i, j) in transitions:
                    ts, tv = 1, 0
                else:
                    ts, tv = 0, 1
                yield i, j, ts, tv


def get_K80_process_definition(kappa):
    exit_rates = np.zeros(4)
    row_states = []
    column_states = []
    unnormalized_rates = []
    for i, j, ts, tv in gen_K80():
        row_states.append([i])
        column_states.append([j])
        rate = kappa * ts + tv
        exit_rates[i] += rate
        unnormalized_rates.append(rate)

    # In the K80 model all exit rates are equal,
    # and each is equal to the expected rate.
    expected_rate = exit_rates.mean()
    assert_allclose(exit_rates, expected_rate)
    transition_rates = [r / expected_rate for r in unnormalized_rates]

    return dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates)


def main():

    # Define the root prior which is uniform for this whole project.
    root_prior = dict(
            states = [[0], [1], [2], [3]],
            probabilities = [0.25, 0.25, 0.25, 0.25])

    # Define the four edges connecting the five nodes of the tree.
    row_nodes = [0, 0, 1, 1]
    column_nodes = [2, 1, 3, 4]
    
    # Define the tree used for simulation in the Minin and Suchard paper.
    # In our case, we will use this to compute log likelihoods for all possible
    # alignment site patterns and to subsequently use these likelihoods
    # as site-specific weights.
    simulation_tree = dict(
            row_nodes = row_nodes,
            column_nodes = column_nodes,
            edge_rate_scaling_factors = [0.3, 0.2, 0.1, 0.1],
            edge_processes = [0, 0, 0, 0])

    # Define the simulation process.
    simulation_process = get_K80_process_definition(4)

    # Define all possible site patterns.
    all_site_patterns = list(itertools.product(range(4), repeat=3))

    # Observation data consisting of all site patterns.
    observed_all_site_patterns = dict(
            nodes = [2, 3, 4],
            variables = [0, 0, 0],
            iid_observations = all_site_patterns)

    # Define the scene associated with the simulation.
    scene = dict(
            node_count = 5,
            process_count = 1,
            state_space_shape = [4],
            tree = simulation_tree,
            root_prior = root_prior,
            process_definitions = [simulation_process],
            observed_data = observed_all_site_patterns)

    # Define the request for per-pattern log likelihoods.
    log_likelihoods_request = {'property' : 'DNNLOGL'}

    # Get the per-pattern log likelihoods.
    j_in = dict(
            scene = scene,
            requests = [log_likelihoods_request])
    j_out = process_json_in(j_in)
    pattern_log_likelihoods = j_out['responses'][0]
    pattern_likelihoods = np.exp(pattern_log_likelihoods)

    print(pattern_likelihoods)
    print(pattern_likelihoods.sum())

    # Get the 

main()
