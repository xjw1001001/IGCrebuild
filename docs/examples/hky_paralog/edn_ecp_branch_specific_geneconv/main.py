"""
This example is based on the tut08 example.

Here we will look at edge-specific expectations.

"""
from __future__ import print_function, division, absolute_import

import functools
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


def gen_geneconv_tau_transition_mask(distn, kappa, tau):
    """
    Yield pairs of multivariate states and the geneconv rate proportion.

    This function depends on the structure and on the spcific parameter values
    of the HKY85+IGC model.

    """
    R, expected_rate = hky(distn, kappa)
    R = R / expected_rate
    for i in range(4):
        for j in range(4):
            if i != j:
                yield [i, j], [i, i], tau / (R[j, i] + tau)
                yield [i, j], [j, j], tau / (R[i, j] + tau)


def gen_heterogeneous_states():
    """
    Yield heterogeneous multivariate states.

    This does not involve edge rate scaling factors
    or process-specific parameter values.
    This function will help compute the proportion of the edge time spent in
    heterogeneous multivariate states.

    """
    for i in range(4):
        for j in range(4):
            if i != j:
                yield [i, j]


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

def get_process_defn_and_prior(distn, kappa, tau):
    triples = list(gen_transitions(distn, kappa, tau))
    rows, cols, transition_rates = zip(*triples)
    process_definition = {
            'row_states' : [list(x) for x in rows],
            'column_states' : [list(x) for x in cols],
            'transition_rates' : list(transition_rates)
            }
    root_prior = {
            "states" : [[0, 0], [1, 1], [2, 2], [3, 3]],
            "probabilities" : list(distn)
            }
    return process_definition, root_prior

def objective_and_gradient(scene, X):
    delta = 1e-8
    distn, kappa, tau, rates, unpacking_cost = unpack(X)
    scene['tree']['edge_rate_scaling_factors'] = rates.tolist()
    log_likelihood_request = {'property' : 'snnlogl'}
    derivatives_request = {'property' : 'sdnderi'}

    # Get the log likelihood and per-edge derivatives.
    # Note that the edge derivatives are of the log likelihood
    # with respect to logs of edge rates, and we will eventually
    # multiply them by -1 to get the gradient of the cost function
    # which we want to minimize rather than the log likelihood function
    # which we want to maximize.
    process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
    scene['root_prior'] = root_prior
    scene['process_definitions'] = [process_defn]
    j_in = {
            'scene' : scene,
            'requests' : [log_likelihood_request, derivatives_request]
            }
    j_out = jsonctmctree.interface.process_json_in(j_in)
    log_likelihood, edge_gradient = j_out['responses']
    cost = -log_likelihood + unpacking_cost

    # For each non-edge-specific parameter get finite-differences
    # approximation of the gradient.
    nedges = len(scene['tree']['row_nodes'])
    nparams = len(X) - nedges
    gradient = []
    for i in range(nparams):
        W = np.copy(X)
        W[i] += delta
        distn, kappa, tau, rates, unpacking_cost = unpack(W)
        process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
        scene['root_prior'] = root_prior
        scene['process_definitions'] = [process_defn]
        j_in = {
                'scene' : scene,
                'requests' : [log_likelihood_request]
                }
        j_out = jsonctmctree.interface.process_json_in(j_in)
        ll = j_out['responses'][0]
        c = -ll + unpacking_cost
        slope = (c - cost) / delta
        gradient.append(slope)
    gradient.extend([-x for x in edge_gradient])
    gradient = np.array(gradient)

    # Return cost and gradient.
    return cost, gradient

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

    # Hard-code the edges of the tree.
    # Note that zip(*x) transposes an array, providing tuples
    # which need to be converted to lists.
    edges = [
            [5, 0], [5, 6],
            [6, 1], [6, 7],
            [7, 2], [7, 8],
            [8, 3], [8, 4]]
    row_nodes, column_nodes = zip(*edges)
    row_nodes = list(row_nodes)
    column_nodes = list(column_nodes)
    node_count = 9
    edge_count = 8
    assert_equal(len(edges), edge_count)
    assert_equal(len(row_nodes), edge_count)
    assert_equal(len(column_nodes), edge_count)

    distn = [0.25, 0.25, 0.25, 0.25]
    edge_rates = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    kappa = 2.0
    tau = 3.0
    process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
    scene = {
            "node_count" : node_count,
            "process_count" : 1,
            "state_space_shape" : [4, 4],
            "tree" : {
                "row_nodes" : row_nodes,
                "column_nodes" : column_nodes,
                "edge_rate_scaling_factors" : edge_rates,
                "edge_processes" : [0, 0, 0, 0, 0, 0, 0, 0]
                },
            "root_prior" : root_prior,
            "process_definition" : process_defn,
            "observed_data" : {
                "nodes" : nodes,
                "variables" : variables,
                "iid_observations" : columns
                }
            }

    X = pack(distn, kappa, tau, edge_rates)
    f = functools.partial(objective_and_gradient, scene)
    result = minimize(f, X, jac=True, method='L-BFGS-B')
    print('final value of objective function:', result.fun)
    distn, kappa, tau, edge_rates, unpacking_cost = unpack(result.x)
    print('nucleotide distribution:')
    for nt, p in zip('ACGT', distn):
        print('  ', nt, ':', p)
    print('kappa:', kappa)
    print('tau:', tau)
    print('edge rate scaling factors:')
    for r in edge_rates:
        print('  ', r)


    # Now for the next phase,
    # compute posterior edge-specific expectations.

    # First, update the scene with the maximum likelihood estimates.
    process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
    scene['process_definitions'] = [process_defn]
    scene['root_prior'] = root_prior
    scene['tree']['edge_specific_scaling_factors'] = edge_rates

    # Prepare the log likelihood request as a control.
    log_likelihood_request = dict(
            property = 'SNNLOGL')

    # Prepare the request for the proportion of time
    # spent in a heterogeneous state on each edge.
    # Do not reduce over states, and do not yet involve
    # the edge-specific rate scaling factors.
    heterogeneous_states = list(gen_heterogeneous_states())
    dwell_request = dict(
            property = 'SDWDWEL',
            state_reduction = dict(
                states = heterogeneous_states,
                weights = [2] * len(heterogeneous_states)))

    # Prepare the request for the gene conversions on edges.
    row_states = []
    column_states = []
    proportions = []
    for info in gen_geneconv_tau_transition_mask(distn, kappa, tau):
        row_state, column_state, proportion = info
        row_states.append(row_state)
        column_states.append(column_state)
        proportions.append(proportion)
    transition_request = dict(
            property = 'SDNTRAN',
            transition_reduction = dict(
                row_states = row_states,
                column_states = column_states,
                weights = proportions))

    # Process the requests.
    j_in = dict(
        scene = scene,
        requests = [
            log_likelihood_request,
            dwell_request,
            transition_request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    responses = j_out['responses']

    ll_response, dwell_response, transition_response = responses

    print('edge specific posterior estimates of tau:')
    for dwel, tran, rate in zip(
            dwell_response, transition_response, edge_rates):
        edge_specific_tau = tran / (dwel * rate)
        print('  ', edge_specific_tau)


main()
