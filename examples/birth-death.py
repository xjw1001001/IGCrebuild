"""
"""
from __future__ import print_function, division

from jsonctmctree.interface import process_json_in

def gen_speciation_triples(lam, n):
    for i in range(1, n-1):
        yield [1, i], [1, i+1], i*lam

def gen_extinction_triples(mu, n):
    yield [1, 1], [0, 0], mu
    for i in range(2, n):
        yield [1, i], [1, i-1], i*mu

def create_j_in(n):

    # speciation rate
    lam = 0.5

    # extinction rate
    mu = 0.3

    # Edge rate scaling factors.
    # In some contexts, these should be treated as times.
    edge_rates = [10, 5, 5, 5]

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
    # The weighted sum over states consists of the 
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

    # This is an upper bound
    # on the probability to have exceeded the population cap.
    misspecification_request = dict(
            property='ssntran',
            transition_reduction = dict(
                row_states = [[1, n-2]],
                column_states = [[1, n-1]],
                weights = [1]))

    j_in = dict(
            scene=scene,
            requests=[
                extinction_request,
                extant_request,
                dwell_request,
                misspecification_request,
                ])

    return j_in

def main():
    j_in = create_j_in(1000)
    print(process_json_in(j_in))

main()
