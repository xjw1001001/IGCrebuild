from __future__ import print_function, division, absolute_import

import json

from jsonctmctree.extras import optimize_em
from jsonctmctree.interface import process_json_in


def main():
    with open('in01.json') as fin:
        j_in = json.load(fin)
    scene = j_in['scene']
    observation_reduction = j_in['requests'][0]['observation_reduction']
    node_count = scene['node_count']
    edge_count = node_count - 1

    # These starting points for EM are OK.
    rates = [0.001 for r in range(edge_count)]
    #rates = [0.01 for r in range(edge_count)]
    #rates = [0.1 for r in range(edge_count)]

    # These starting points are questionable.
    #rates = [0.15 for r in range(edge_count)]
    #rates = [0.2 for r in range(edge_count)]

    # These starting points are not really feasible.
    #rates = [0.5 for r in range(edge_count)]
    #rates = [0.95 for r in range(edge_count)]
    #rates = [1 for r in range(edge_count)]

    # Initialize rates.
    scene['tree']['edge_rate_scaling_factors'] = rates

    # Update rates according to EM.
    rates = optimize_em(j_in['scene'], observation_reduction, 3)

    # Show the log likelihood
    scene['tree']['edge_rate_scaling_factors'] = rates
    j_in = dict(
            scene = scene,
            requests = [dict(
                property = 'WNNLOGL',
                observation_reduction = observation_reduction)])
    ll = process_json_in(j_in)['responses'][0]
    print(ll)

if __name__ == '__main__':
    main()
