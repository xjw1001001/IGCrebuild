"""
"""
from __future__ import print_function, division

from math import exp, expm1
import argparse

from jsonctmctree.interface import process_json_in

# high population boundary behaviors
WRAP = 'wrap'
ABSORB = 'absorb'
BLOCK = 'block'

def gen_speciation_triples(lam, n, boundary):
    for i in range(1, n-1):
        yield [1, i], [1, i+1], i*lam
    if boundary == WRAP:
        yield [1, n-1], [0, 0], (n-1)*lam

def gen_extinction_triples(mu, n, boundary):
    yield [1, 1], [0, 0], mu
    for i in range(2, n-1):
        yield [1, i], [1, i-1], i*mu
    if boundary in {WRAP, BLOCK}:
        yield [1, n-1], [1, n-2], (n-1)*mu

def get_scene(edge_rates, mu, lam, n, boundary):
    leaves = [2, 3, 4]
    tree = dict(
            row_nodes = [0, 0, 1, 1],
            column_nodes = [1, 2, 3, 4],
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0, 0, 0, 0])
    speciation_triples = list(gen_speciation_triples(lam, n, boundary))
    extinction_triples = list(gen_extinction_triples(mu, n, boundary))
    triples = speciation_triples + extinction_triples
    row_states, column_states, transition_rates = zip(*triples)
    scene = dict(
            node_count = 5,
            process_count = 1,
            state_space_shape = [2, n],
            tree = tree,
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
    return scene

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=100, help='population size')
    parser.add_argument('--mu', type=float, default=0.3, help='extinction')
    parser.add_argument('--lam', type=float, default=0.5, help='speciation')
    args = parser.parse_args()

    n = args.n
    mu = args.mu
    lam = args.lam

    edge_rates = [5, 10, 5, 5]
    wrap_scene = get_scene(edge_rates, mu, lam, n, WRAP)
    absorb_scene = get_scene(edge_rates, mu, lam, n, ABSORB)

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

    # Compute only the likelihood for the absorbing high population boundary.
    j_out = process_json_in(dict(scene=absorb_scene, requests=[logl_request]))
    absorb_likelihood = exp(j_out['responses'][0])

    # Compute more stuff for the wrapping boundary.
    j_in = dict(
            scene=wrap_scene,
            requests=[
                logl_request,
                extinction_request,
                extant_request,
                dwell_request,
                ])
    j_out = process_json_in(j_in)

    logl, extinction, extant, dwell = j_out['responses']
    wrap_likelihood = exp(logl)
    print('likelihood:', wrap_likelihood)
    print('upper bound likelihood for unbounded population:', absorb_likelihood)
    print('expected number of extinctions:', extinction)
    print('expected number of extant lineages at each node:')
    for i, x in enumerate(extant):
        print(i, ':', x)
    print('expected total size of the gene tree:', dwell)

main()
