"""
The NCD example from marca_doc.ps.

"""
from __future__ import print_function, division, absolute_import

import numpy as np

import jsonctmctree.interface


def gen_discrete_simplex(n, k):
    """
    n : the population size
    k : the number of variables
    yield suffixes
    """
    if not n:
        yield [0] * k
    elif k == 1:
        yield [n]
    else:
        for i in range(n+1):
            prefix = [i]
            for suffix in gen_discrete_simplex(n - i, k - 1):
                yield prefix + suffix


def gen_triples(n):
    # initial state, final state, rate
    states = list(gen_discrete_simplex(n, 4))
    for sa in states:
        CPU, SM, FD, TERM = sa
        if CPU:
            SM_rate = 100 * ((CPU + SM + FD) / 128)**1.5
            if SM < 10:
                yield sa, [CPU-1, SM+1, FD, TERM], SM_rate
            if FD < 10:
                yield sa, [CPU-1, SM, FD+1, TERM], 0.05
            if TERM < 10:
                yield sa, [CPU-1, SM, FD, TERM+1], 0.002
        if CPU < 10:
            if SM:
                yield sa, [CPU+1, SM-1, FD, TERM], 0.2
            if FD:
                yield sa, [CPU+1, SM, FD-1, TERM], 1/30
            if TERM:
                yield sa, [CPU+1, SM, FD, TERM-1], 0.0001 * TERM


def get_process_definition():
    triples = list(gen_triples(10))
    row_states, column_states, rates = zip(*triples)
    return dict(
            row_states = list(row_states),
            column_states = list(column_states),
            transition_rates = list(rates))


def get_scene(t, observations):
    # t is the time.

    tree = dict(
            row_nodes = [0],
            column_nodes = [1],
            edge_rate_scaling_factors = [t],
            edge_processes = [0])

    root_prior = dict(
            states = [
                [4, 0, 0, 6],
                [4, 1, 0, 5],
                [4, 0, 1, 5],
                [4, 2, 0, 4]],
            probabilities = [
                0.25, 0.25, 0.25, 0.25])

    # The first three observations are specifically requested.
    # The sum of the last four observations gives the total
    observed_data = dict(
            nodes = [1, 1, 1, 1],
            variables = [0, 1, 2, 3],
            iid_observations = observations)

    process_definition = get_process_definition()

    # Assemble the scene.
    scene = dict(
            node_count = 2,
            process_count = 1,
            state_space_shape = [11, 11, 11, 11],
            tree = tree,
            root_prior = root_prior,
            process_definitions = [process_definition],
            observed_data = observed_data)
    
    return scene


def main():

    # Define some requests.
    # These include the log likelihood
    # and the sum of transition count expectations.
    log_likelihood_request = dict(property = 'DNNLOGL')

    # Just get likelihoods for all possible observations.
    observations = list(gen_discrete_simplex(10, 4))
    obs_to_idx = dict((tuple(x), i) for i, x in enumerate(observations))

    for t in 10, 25:
        print('t:', t)
        print('-----')
        print()

        # Get the per-state log likelihoods.
        scene = get_scene(t, observations)
        j_in = {
                "scene" : scene,
                "requests" : [log_likelihood_request],
                }
        j_out = jsonctmctree.interface.process_json_in(j_in)
        log_likelihoods = j_out['responses'][0]
        likelihoods = np.exp(log_likelihoods)

        print('custom state probability requests:')
        for obs in (
                (0, 0, 0, 10),
                (6, 0, 0, 4),
                (4, 0, 0, 6)):
            p = likelihoods[obs_to_idx[obs]]
            print(list(obs), p, sep='\t')
        print()

        print('cumulative probabilities for constrained states:')
        cumulative = 0
        for obs in (
                (4, 0, 0, 6),
                (4, 1, 0, 5),
                (4, 0, 1, 5),
                (4, 2, 0, 4),
                (4, 1, 1, 4),
                (4, 0, 2, 4),
                (4, 3, 0, 3),
                (4, 2, 1, 3),
                (4, 1, 2, 3),
                (4, 0, 3, 3)):
            p = likelihoods[obs_to_idx[obs]]
            cumulative += p
            print(list(obs), p, cumulative, sep='\t')
        print()
        
        for variable in range(4):
            print('marginal distribution of variable %d:' % (variable + 1))
            distn = np.zeros(11)
            for obs in observations:
                p = likelihoods[obs_to_idx[tuple(obs)]]
                distn[obs[variable]] += p
            for i, m in enumerate(distn):
                print(i, m, sep='\t')
            mu = np.dot(distn, range(11))
            var = np.dot(distn, np.square(np.arange(11) - mu))
            std = np.sqrt(var)
            print('expected value:', mu)
            print('standard deviation:', std)
            print()
        print()


main()
