import json
import copy

import jsonctmctree.interface


def main():
    n = 20
    with open('in01.json') as fin:
        j_in = json.load(fin)
    j_out = jsonctmctree.interface.process_json_in(j_in)
    print(j_out)

    # do a few EM iterations
    nsites = sum(j_in['requests'][0]['observation_reduction']['weights'])
    for i in range(10):

        # update branch lengths
        expected_rate = n - 1
        expected_events_per_edge = j_out['responses'][1]
        j_in['scene']['tree']['edge_rate_scaling_factors'] = [
                x / (expected_rate * nsites) for x in expected_events_per_edge]
        j_out = jsonctmctree.interface.process_json_in(j_in)

        print(j_out)

if __name__ == '__main__':
    main()
