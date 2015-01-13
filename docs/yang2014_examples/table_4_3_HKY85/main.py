# NOTE: the iterative algorithm used by this script may not be a true EM,
#       and it may not converge to the maximum likelihood estimates.

from __future__ import print_function, division, absolute_import

from functools import partial
from collections import defaultdict
import copy
import json

import pyparsing

import numpy as np
from numpy.testing import assert_equal
import scipy.optimize

import jsonctmctree.interface

s_tree = """(kiwi_fruit, ((((agave, garlic), rice), black_pepper),
((cabbage, cotton), (cucumber, walnut))), (sunflower, (tomato, tobacco)))"""


# An observed state of -1 means completely missing data.
_nt_to_state = {
        'A' : [1, 0, 0, 0],
        'C' : [0, 1, 0, 0],
        'G' : [0, 0, 1, 0],
        'T' : [0, 0, 0, 1],
        '?' : [-1, -1, -1, -1],
        '-' : [-1, -1, -1, -1],
        'N' : [-1, -1, -1, -1],
        'M' : [-1, -1, 0, 0],
        'R' : [-1, 0, -1, 0],
        'W' : [-1, 0, 0, -1],
        'S' : [0, -1, -1, 0],
        'Y' : [0, -1, 0, -1],
        'K' : [0, 0, -1, -1],
        }

def _help_build_tree(parent, root, node, name_to_node, edges):
    if parent is not None:
        edges.append((parent, node))
    neo = node + 1
    if isinstance(root, basestring):
        name_to_node[root] = node
    else:
        for element in root:
            neo = _help_build_tree(node, element, neo, name_to_node, edges)
    return neo

def get_tree_info():
    # Return a dictionary mapping name to node index,
    # and return a list of edges as ordered pairs of node indices.
    tree = s_tree.replace(',', ' ')
    nestedItems = pyparsing.nestedExpr(opener='(', closer=')')
    tree = (nestedItems + pyparsing.stringEnd).parseString(tree).asList()[0]
    name_to_node = {}
    edges = []
    _help_build_tree(None, tree, 0, name_to_node, edges)
    return name_to_node, edges


def get_maximum_likelihood_pi(character_to_count):
    # Use EM to get the maximum likelihood nucleotide frequency distribution.
    # This gives a distribution similar to the one computed by PAML.
    informative_state_to_indicators = dict()
    for nt, arr in _nt_to_state.items():
        indicators = np.absolute(arr)
        if indicators.sum() < 4:
            informative_state_to_indicators[nt] = indicators
    # Use a few iterations.
    # Begin with a uniform distribution for the first iteration.
    pi = np.ones(4) / 4
    for i in range(10):
        weights = np.zeros(4)
        for character, count in character_to_count.items():
            indicators = informative_state_to_indicators.get(character, None)
            if indicators is not None:
                p = pi * indicators
                weights = weights + count * (p / p.sum())
        pi = weights / weights.sum()
    # Return the EM estimate of the maximum likelihood
    # probability distribution over nucleotides.
    return pi


def get_nucleotide_alignment_info(name_to_node):
    # Return an ordered list of observable node indices,
    # and return the array of iid observations.
    # Some of the states are more informative than others.
    # Eventually use this information to estimate nuceotide frequencies
    # with maximum likelihood by using EM to deal with the
    # partially informative sites.
    #nt_to_count_vector = {}
    #for k, state in _nt_to_state.items():
        #v = np.absolute(state)
        #nt_to_count_vector[k] = v / v.sum()
    nt_to_idx = {c : i for i, c in enumerate('ACGT')}
    sequences = []
    nodes = []
    variables = []
    acgt_counts = np.zeros(4)
    character_to_count = defaultdict(int)
    # Get counts of each completely or partially informative nucleotide.
    with open('vegetables.rbcL.txt') as fin:
        lines = fin.readlines()
        header = lines[0]
        for line in lines[1:]:
            name, sequence = line.strip().split()
            for c in sequence:
                character_to_count[c] += 1
                #acgt_counts += nt_to_count_vector[c]
                idx = nt_to_idx.get(c, None)
                if idx is not None:
                    acgt_counts[idx] += 1
            node = name_to_node[name]
            nodes.extend([node]*4)
            sequences.append([_nt_to_state[c] for c in sequence])
            variables.extend([0, 1, 2, 3])
    # arr[node, site, variable]
    arr = np.array(sequences, dtype=int)
    arr = np.transpose(arr, (1, 0, 2))
    iid_observations = np.reshape(arr, (arr.shape[0], -1))
    iid_observations = iid_observations.tolist()

    # Compute pi empirically.
    #pi = (acgt_counts / acgt_counts.sum()).tolist()
    # Or use the values from the textbook.
    #pi = [0.2754, 0.1927, 0.2452, 0.2867]
    # Or get maximum likelihood estimates with EM in the manner of PAML.
    pi = get_maximum_likelihood_pi(character_to_count)

    return nodes, variables, iid_observations, pi


def get_expected_rate(pi, kappa):
    raw_exit_rates = get_unnormalized_exit_rates(pi, kappa)
    return np.dot(pi, raw_exit_rates)


def get_exit_rates(pi, kappa):
    raw_exit_rates = get_unnormalized_exit_rates(pi, kappa)
    expectation = np.dot(pi, raw_exit_rates)
    return [r / expectation for r in raw_exit_rates]


def get_unnormalized_exit_rates(pi, kappa):
    exit_rates = [0, 0, 0, 0]
    for i, j, ts, tv in gen_hky():
        rate = (kappa * ts + tv) * pi[j]
        exit_rates[i] += rate
    return exit_rates

def get_unnormalized_ts_exits(pi, kappa):
    exit_rates = [0, 0, 0, 0]
    for i, j, ts, tv in gen_hky():
        if ts:
            rate = (kappa * ts + tv) * pi[j]
            exit_rates[i] += rate
    return exit_rates

def get_unnormalized_tv_exits(pi, kappa):
    exit_rates = [0, 0, 0, 0]
    for i, j, ts, tv in gen_hky():
        if tv:
            rate = (kappa * ts + tv) * pi[j]
            exit_rates[i] += rate
    return exit_rates


def get_process_definition(pi, kappa):

    # Put together the multivariate rate matrix,
    # normalized by dividing by the expected rate.
    expected_rate = get_expected_rate(pi, kappa)
    ident = np.identity(4, dtype=int).tolist()
    row_states = []
    column_states = []
    rates = []
    for i, j, ts, tv in gen_hky():
        rate = (kappa * ts + tv) * pi[j]
        row_states.append(ident[i])
        column_states.append(ident[j])
        rates.append(rate / expected_rate)

    # Assemble and return the process definition.
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = rates)
    return process_definition


def gen_hky():
    # Yield a tuple for each state transition.
    # Currently, this tuple consists of the univariate initial state,
    # the univariate final state, a ts indicator, and a tv indicator.
    # ts: A<->G, C<->T
    ts_pairs = ((0, 2), (2, 0), (1, 3), (3, 1))
    for i in range(4):
        for j in range(4):
            if i != j:
                ts = 1 if (i, j) in ts_pairs else 0
                tv = 1 - ts
                yield i, j, ts, tv


def get_requests(edge_rates, pi, kappa):

    # Precompute the structure of the sparse matrix.
    edge_count = len(edge_rates)
    ident = np.identity(4, dtype=int).tolist()
    hky = list(gen_hky())

    # Define the exit rates given the current parameter values.
    # This includes normalization by dividing by the expected rate.
    exit_rates = get_exit_rates(pi, kappa)

    # Define the log likelihood request.
    log_likelihood_request = {"property" : "SNNLOGL"}

    # Define the requests for expecatations that are used
    # to update the branch length parameter estimates.
    per_edge_opportunity_request = dict(
            property = "SDWDWEL",
            state_reduction = dict(
                states = ident,
                weights = exit_rates))
    per_edge_change_request = dict(
            property = "SDNTRAN",
            transition_reduction = dict(
                row_states = [ident[i] for i, j, ts, tv in hky],
                column_states = [ident[j] for i, j, ts, tv in hky],
                weights = [1] * len(hky)))

    # Define the requests for expectations that are used
    # to update the kappa estimates.

    ts_row_states = [ident[i] for i, j, ts, tv in hky if ts]
    ts_column_states = [ident[j] for i, j, ts, tv in hky if ts]
    ts_rates = [kappa * pi[j] for i, j, ts, tv in hky if ts]

    tv_row_states = [ident[i] for i, j, ts, tv in hky if tv]
    tv_column_states = [ident[j] for i, j, ts, tv in hky if tv]
    tv_rates = [pi[j] for i, j, ts, tv in hky if tv]


    # Get unnormalized ts and tv exit rates.
    ts_exits = get_unnormalized_ts_exits(pi, kappa)
    tv_exits = get_unnormalized_tv_exits(pi, kappa)
    edge_reduction = dict(
            edges = range(edge_count),
            weights = edge_rates)

    ts_opportunity_request = dict(
            property = "SWWDWEL",
            edge_reduction = edge_reduction,
            state_reduction = dict(
                states = ident,
                weights = ts_exits))
    tv_opportunity_request = dict(
            property = "SWWDWEL",
            edge_reduction = edge_reduction,
            state_reduction = dict(
                states = ident,
                weights = tv_exits))
    ts_change_request = dict(
            property = "SWNTRAN",
            edge_reduction = edge_reduction,
            transition_reduction = dict(
                row_states = ts_row_states,
                column_states = ts_column_states,
                weights = ts_rates))
    tv_change_request = dict(
            property = "SWNTRAN",
            edge_reduction = edge_reduction,
            transition_reduction = dict(
                row_states = tv_row_states,
                column_states = tv_column_states,
                weights = tv_rates))

    return [
        log_likelihood_request,
        per_edge_opportunity_request,
        per_edge_change_request,
        ts_opportunity_request,
        tv_opportunity_request,
        ts_change_request,
        tv_change_request]

def pack(edge_rates, kappa):
    return np.log(list(edge_rates) + [kappa])

def unpack(X):
    Y = np.exp(X)
    edge_rates = Y[:-1].tolist()
    kappa = Y[-1]
    return edge_rates, kappa

def objective(scene, pi, X):
    scene = copy.deepcopy(scene)
    edge_rates, kappa = unpack(X)
    scene['tree']['edge_rate_scaling_factors'] = edge_rates
    scene['process_definitions'] = [get_process_definition(pi, kappa)]
    j_in = dict(
            scene = scene,
            requests = [{"property" : "SNNLOGL"}])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    log_likelihood = j_out['responses'][0]
    return -log_likelihood


def main():

    # Read the tree.
    name_to_node, edges = get_tree_info()
    edge_count = len(edges)
    node_count = edge_count + 1

    # Read the alignment.
    info = get_nucleotide_alignment_info(name_to_node)
    nodes, variables, iid_observations, empirical_pi = info
    nsites = len(iid_observations)

    # Initialize some guesses.
    edge_rates = [0.01] * edge_count
    pi = empirical_pi
    kappa = 2.0
    #kappa = 3.620

    # Define the tree component of the scene
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0] * edge_count)

    # Define the structure of the hky model.
    hky = list(gen_hky())

    # Define the root distribution.
    ident = np.identity(4, dtype=int).tolist()
    root_prior = dict(
            states = ident,
            probabilities = pi)

    # Define the observed data.
    observed_data = dict(
            nodes = nodes,
            variables = variables,
            iid_observations = iid_observations)

    # Assemble the scene.
    scene = dict(
            node_count = node_count,
            process_count = 1,
            state_space_shape = [2, 2, 2, 2],
            tree = tree,
            root_prior = root_prior,
            observed_data = observed_data)

    arr = []
    j_out = None
    for i in range(10):

        # if j_out is available, recompute kappa and edge rates
        if j_out is not None:
            responses = j_out['responses']
            (
                    ll,
                    per_edge_opportunity,
                    per_edge_change,
                    ts_opportunity,
                    tv_opportunity,
                    ts_change,
                    tv_change) = responses
            edge_rates = []
            for change, dwell in zip(per_edge_change, per_edge_opportunity):
                edge_rates.append(change / dwell)
            kappa = (ts_change / ts_opportunity) / (tv_change / tv_opportunity)

        process_definition = get_process_definition(pi, kappa)
        j_in = dict(scene = scene)
        j_in['scene']['tree']['edge_rate_scaling_factors'] = edge_rates
        j_in['requests'] = get_requests(edge_rates, pi, kappa)
        j_in['scene']['process_definitions'] = [process_definition]
        j_out = jsonctmctree.interface.process_json_in(j_in)
        arr.append(copy.deepcopy(j_out))

        print(j_out)
        print(kappa)

    # Improve the estimates using a numerical search.
    f = partial(objective, scene, pi)
    x0 = pack(edge_rates, kappa)
    result = scipy.optimize.minimize(f, x0, method='L-BFGS-B')
    xopt = result.x
    print(result)
    edge_rates, kappa = unpack(xopt)
    print(edge_rates)
    print(kappa)

    #print(json.dumps(arr, indent=4))


main()
