from __future__ import print_function, division, absolute_import

import json
import copy

import jsonctmctree.interface

def main():
    n = 20
    expected_rate = n - 1
    j_out = None
    with open('in01.json') as fin:
        j_in = json.load(fin)
    nsites = sum(j_in['requests'][0]['observation_reduction']['weights'])
    arr = []
    for i in range(4):
        if j_out is None:
            j_out = jsonctmctree.interface.process_json_in(j_in)
        else:
            expected_events = j_out['responses'][1]
            j_in['scene']['tree']['edge_rate_scaling_factors'] = [
                    x / (expected_rate * nsites) for x in expected_events]
            j_out = jsonctmctree.interface.process_json_in(j_in)
        arr.append(copy.deepcopy(j_out))
    print(json.dumps(arr))

if __name__ == '__main__':
    main()
