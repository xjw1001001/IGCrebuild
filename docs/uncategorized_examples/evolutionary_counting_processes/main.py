"""
This is really a K80 model because the nucleotide distribution is uniform.

Reproduce an example from Minin and Suchard
"Counting labeled transitions in continuous-time Markov models of evolution."

"""
from __future__ import print_function, division

import itertools
import copy

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


def run_analysis(scene, pattern_likelihoods, kappa):

    # Copy the scene because we are going to do some surgery.
    scene = copy.deepcopy(scene)

    # Change the edge rate scaling factors.
    scene['tree']['edge_rate_scaling_factors'] = [0.28, 0.21, 0.12, 0.09]

    # Define the model to be used for the analysis.
    analysis_process = get_K80_process_definition(kappa)
    scene['process_definitions'] = [analysis_process]

    # Define the observation reduction to be used in both requests.
    npatterns = len(pattern_likelihoods)
    observation_reduction = dict(
            observation_indices = range(npatterns),
            weights = pattern_likelihoods)

    # Request nucleotide transition count expectations.
    ts_transition_request = dict(
            property = 'WSNTRAN',
            observation_reduction = observation_reduction,
            transition_reduction = dict(
                row_states = [[i] for i, j, ts, tv in gen_K80() if ts],
                column_states = [[j] for i, j, ts, tv in gen_K80() if ts],
                weights = [1 for i, j, ts, tv in gen_K80() if ts]))

    # Request nucleotide transversion count expectations.
    tv_transition_request = dict(
            property = 'WSNTRAN',
            observation_reduction = observation_reduction,
            transition_reduction = dict(
                row_states = [[i] for i, j, ts, tv in gen_K80() if tv],
                column_states = [[j] for i, j, ts, tv in gen_K80() if tv],
                weights = [1 for i, j, ts, tv in gen_K80() if tv]))

    # Define the requests.
    requests = [ts_transition_request, tv_transition_request]

    # Run the analysis.
    j_in = dict(
            scene = scene,
            requests = requests)
    j_out = process_json_in(j_in)
    print(j_out)



def main():

    # The root prior for every process in this project is uniform.
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

    # The sum of likelihoods over all patterns should be 1.
    assert_allclose(pattern_likelihoods.sum(), 1)

    # Run an analysis for each of a few values of kappa.
    for kappa in 1, 2, 4:
        run_analysis(scene, pattern_likelihoods.tolist(), kappa)


main()
