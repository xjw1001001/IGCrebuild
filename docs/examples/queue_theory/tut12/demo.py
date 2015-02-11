from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt

from jsonctmctree.interface import process_json_in

def main():
    ts = np.linspace(1e-5, 30, 100)
    n = len(ts)
    j_in = {
        "scene" : {
            "node_count" : n+1,
            "process_count" : 1,
            "state_space_shape" : [4],
            "tree" : {
                "row_nodes" : [n]*n,
                "column_nodes" : range(n),
                "edge_rate_scaling_factors" : (0.5 * ts).tolist(),
                "edge_processes" : [0]*n
            },
            "root_prior" : {
                "states" : [[0]],
                "probabilities" : [1]
            },
            "process_definitions" : [{
                "row_states" : [[0], [1], [2]],
                "column_states" : [[1], [2], [3]],
                "transition_rates" : [1, 2, 3]
            }],
            "observed_data" : {
                "nodes" : [],
                "variables" : [],
                "iid_observations" : [[]]
            }
        },
        "requests" : [{"property" : "SDDDWEL"}]
    }
    j_out = process_json_in(j_in)
    a, b, c, d = zip(*j_out['responses'][0])
    lines = plt.plot(
            ts, a, 'blue',
            ts, b, 'green',
            ts, c, 'red',
            ts, d, 'skyblue')
    plt.ylabel("Time-averaged Expected sojourn time")
    plt.xlabel("Time")
    plt.legend(
            lines,
            ('State 1', 'State 2', 'State 3', 'State 4 (absorbing)'),
            loc='center right')
    plt.savefig('out00.svg', transparent=True)

main()
