from __future__ import print_function, division, absolute_import

import functools
import itertools
import copy
import json

import numpy as np
from numpy.testing import assert_equal
from scipy.misc import logsumexp
from scipy.optimize import minimize

import jsonctmctree.interface

def hky(distn, k):
    R = np.array([
        [0, 1, k, 1],
        [1, 0, 1, k],
        [k, 1, 0, 1],
        [1, k, 1, 0],
        ]) * distn
    return R, R.sum(axis=1).dot(distn)

def gen_transitions(distn, kappa, tau):
    R, expected_rate = hky(distn, kappa)
    R = R / expected_rate
    for i in range(4):
        for j in range(4):
            if i == j:
                for k in range(4):
                    if i != k:
                        yield (i, j), (k, j), R[i, k]
                    if j != k:
                        yield (i, j), (i, k), R[j, k]
            else:
                yield (i, j), (i, i), R[j, i] + tau
                yield (i, j), (j, j), R[i, j] + tau
                for k in range(4):
                    if i != k and j != k:
                        yield (i, j), (k, j), R[i, k]
                        yield (i, j), (i, k), R[j, k]

def pack(distn, kappa, tau, rates):
    return np.log(np.concatenate([distn, [kappa, tau], rates]))

def unpack(X):
    lse = logsumexp(X[0:4])
    unpacking_cost = lse * lse
    distn = np.exp(X[0:4] - lse)
    kappa, tau = np.exp(X[4:6])
    rates = np.exp(X[6:])
    return distn, kappa, tau, rates, unpacking_cost

def objective(scene, X):
    distn, kappa, tau, rates, unpacking_cost = unpack(X)
    scene['root_prior']['probabilities'] = distn.tolist()
    scene['tree']['edge_rate_scaling_factors'] = rates.tolist()
    triples = list(gen_transitions(distn, kappa, tau))
    rows, cols, transition_rates = zip(*triples)
    process_definition = {
            'row_states' : [list(x) for x in rows],
            'column_states' : [list(x) for x in cols],
            'transition_rates' : list(transition_rates)
            }
    scene['process_definitions'] = [process_definition]
    request = {'property' : 'snnlogl'}
    j_in = {'scene' : scene, 'requests' : [request]}
    j_out = jsonctmctree.interface.process_json_in(j_in)
    log_likelihood = j_out['responses'][0]
    cost = -log_likelihood + unpacking_cost
    return cost

def main():
    name_to_node = {
            'tamarin' : 0,
            'macaque' : 1,
            'orangutan' : 2,
            'chimpanzee' : 3,
            'gorilla' : 4}
    paralog_to_variable = {
            'ecp' : 0,
            'edn' : 1}
    nodes = []
    variables = []
    rows = []
    with open('paralogs.fasta') as fin:
        while True:
            line = fin.readline().strip().lower()
            if not line:
                break
            name = line[1:-3]
            paralog = line[-3:]
            seq = fin.readline().strip()
            row = ['ACGT'.index(x) for x in seq]
            nodes.append(name_to_node[name])
            variables.append(paralog_to_variable[paralog])
            rows.append(row)
    columns = [list(x) for x in zip(*rows)]

    # This is a partial scene, missing the root distribution,
    # the process definition, and the observed data.
    distn = [0.25, 0.25, 0.25, 0.25]
    rates = [1, 1, 1, 1, 1, 1, 1, 1]
    scene = {
            "node_count" : 9,
            "process_count" : 1,
            "state_space_shape" : [4, 4],
            "tree" : {
                "row_nodes" : [5, 5, 6, 6, 7, 7, 8, 8],
                "column_nodes" : [0, 6, 1, 7, 2, 8, 3, 4],
                "edge_rate_scaling_factors" : rates,
                "edge_processes" : [0, 0, 0, 0, 0, 0, 0, 0]
                },
            "root_prior" : {
                "states" : [[0, 0], [1, 1], [2, 2], [3, 3]],
                "probabilities" : distn
                },
            "observed_data" : {
                "nodes" : nodes,
                "variables" : variables,
                "iid_observations" : columns
                }
            }

    X = pack(distn, 2.0, 3.0, rates)
    f = functools.partial(objective, scene)
    result = minimize(f, X, method='L-BFGS-B')
    print('final value of objective function:', result.fun)
    distn, kappa, tau, rates, unpacking_cost = unpack(result.x)
    print('nucleotide distribution:')
    for nt, p in zip('ACGT', distn):
        print('  ', nt, ':', p)
    print('kappa:', kappa)
    print('tau:', tau)
    print('edge rate scaling factors:')
    for r in rates:
        print('  ', r)

main()

