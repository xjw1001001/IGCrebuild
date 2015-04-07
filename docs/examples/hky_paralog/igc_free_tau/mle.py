from __future__ import print_function, division, absolute_import

import functools
import itertools
import copy
import json

import numpy as np
from numpy.testing import assert_equal
from scipy.misc import logsumexp
from scipy.optimize import minimize
from scipy.special import logit, expit

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


def pack_acgt(pi):
    a, c, g, t = pi
    ag = a+g  # purines
    ct = c+t  # pyrimidines
    a_div_ag = a / ag
    c_div_ct = c / ct
    return logit([ag, a_div_ag, c_div_ct])


def unpack_acgt(packed_acgt):
    ag, a_div_ag, c_div_ct = expit(packed_acgt)
    ct = 1 - ag
    a = a_div_ag * ag
    g = ag - a
    c = c_div_ct * ct
    t = ct - c
    return np.array([a, c, g, t])


def pack_global_params(pi, kappa, tau):
    return np.concatenate([
        pack_acgt(pi),
        np.log([kappa, tau])])


def unpack_global_params(X):
    pi = unpack_acgt(X[:3])
    kappa, tau = np.exp(X[3:])
    return pi, kappa, tau


def pack(distn, kappa, tau, rates):
    return np.concatenate((
        pack_global_params(distn, kappa, tau),
        np.log(rates)))


def unpack(X):
    distn, kappa, tau = unpack_global_params(X[:5])
    rates = np.exp(X[5:])
    return distn, kappa, tau, rates


def objective(scene, X):
    distn, kappa, tau, rates = unpack(X)
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
    cost = -log_likelihood
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

    print('number of sites in the alignment:', len(columns))
    print('number of sequences:', len(nodes))

    # Compute the empirical distribution of the nucleotides.
    counts = np.zeros(4)
    for k in np.ravel(columns):
        counts[k] += 1
    empirical_pi = counts / counts.sum()

    distn = empirical_pi
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
    distn, kappa, tau, rates = unpack(result.x)
    print('nucleotide distribution:')
    for nt, p in zip('ACGT', distn):
        print('  ', nt, ':', p)
    print('kappa:', kappa)
    print('tau:', tau)
    print('edge rate scaling factors:')
    for r in rates:
        print('  ', r)

main()

