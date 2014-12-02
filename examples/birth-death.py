"""
"""
from __future__ import print_function, division

from math import exp, expm1
import argparse

from jsonctmctree.interface import process_json_in

def gen_speciation_triples(lam, n):
    for i in range(1, n-1):
        yield [1, i], [1, i+1], i*lam
    yield [1, n-1], [0, 0], (n-1)*lam

def gen_extinction_triples(mu, n):
    yield [1, 1], [0, 0], mu
    for i in range(2, n):
        yield [1, i], [1, i-1], i*mu

def create_j_in(mu, lam, n):

    # Edge rate scaling factors.
    # In some contexts, these should be treated as times.
    edge_rates = [5, 10, 5, 5]

    # Define the process.
    speciation_triples = list(gen_speciation_triples(lam, n))
    extinction_triples = list(gen_extinction_triples(mu, n))
    triples = speciation_triples + extinction_triples
    row_states, column_states, transition_rates = zip(*triples)

    leaves = [2, 3, 4]

    scene = dict(
            node_count = 5,
            process_count = 1,
            state_space_shape = [2, n],
            tree = dict(
                    row_nodes = [0, 0, 1, 1],
                    column_nodes = [1, 2, 3, 4],
                    edge_rate_scaling_factors = edge_rates,
                    edge_processes = [0, 0, 0, 0],
                    ),
            root_prior = dict(
                states = [[1, 1]],
                probabilities = [1.0],
                ),
            process_definitions = [dict(
                row_states = row_states,
                column_states = column_states,
                transition_rates = transition_rates,
                )],
            observed_data = dict(
                nodes = leaves,
                variables = [0, 0, 0],
                iid_observations = [[1, 1, 1]],
                )
            )

    # Log likelihood.
    logl_request = dict(property='snnlogl')

    # Unweighted sum over observations and over edges,
    # and weighted sum over transitions consisting of the unweighted sum
    # over transitions corresponding to extinction events.
    extinction_request = dict(
        property='ssntran',
        transition_reduction = dict(
            row_states = [[1, i] for i in range(2, n)],
            column_states = [[1, i-1] for i in range(2, n)],
            weights = [1]*(n-2),
        ))

    # Unweighted sum over observations, and weighted sum over states.
    # TODO implement node reduction in jsonctmctree?
    extant_request = dict(
        property='snwnode',
        state_reduction = dict(
            states = [[1, i] for i in range(n)],
            weights = range(n),
            ))

    # Unweighted sum over observations, weighted sum over edges,
    # and weighted sum over states.
    dwell_request = dict(
        property='swwdwel',
        edge_reduction = dict(
            edges = [0, 1, 2, 3],
            weights = edge_rates),
        state_reduction = dict(
            states = [[1, i] for i in range(n)],
            weights = range(n)))

    j_in = dict(
            scene=scene,
            requests=[
                logl_request,
                extinction_request,
                extant_request,
                dwell_request,
                ])

    return j_in

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=100, help='population size')
    parser.add_argument('--mu', type=float, default=0.3, help='extinction')
    parser.add_argument('--lam', type=float, default=0.5, help='speciation')
    args = parser.parse_args()
    j_in = create_j_in(args.mu, args.lam, args.n)
    j_out = process_json_in(j_in)
    logl, extinction, extant, dwell = j_out['responses']
    print('likelihood:', exp(logl))
    print('expected number of extinctions:', extinction)
    print('expected number of extant lineages at each node:')
    for i, x in enumerate(extant):
        print(i, ':', x)
    print('expected total size of the gene tree:', dwell)

main()
