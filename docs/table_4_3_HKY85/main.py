from __future__ import print_function, division, absolute_import

import copy
import json

import pyparsing

import numpy as np
from numpy.testing import assert_equal

import jsonctmctree.interface

s_tree = """(kiwi_fruit, ((((agave, garlic), rice), black_pepper),
((cabbage, cotton), (cucumber, walnut))), (sunflower, (tomato, tobacco)))"""

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


def get_nucleotide_alignment_info(name_to_node):
    # Return an ordered list of observable node indices,
    # and return the array of iid observations.
    # An observed state of -1 means completely missing data.
    nt_to_state = {
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
    #nt_to_count_vector = {}
    #for k, state in nt_to_state.items():
        #v = np.absolute(state)
        #nt_to_count_vector[k] = v / v.sum()
    nt_to_idx = {c : i for i, c in enumerate('ACGT')}
    sequences = []
    nodes = []
    variables = []
    acgt_counts = np.zeros(4)
    with open('vegetables.rbcL.txt') as fin:
        lines = fin.readlines()
        header = lines[0]
        for line in lines[1:]:
            name, sequence = line.strip().split()
            for c in sequence:
                #acgt_counts += nt_to_count_vector[c]
                idx = nt_to_idx.get(c, None)
                if idx is not None:
                    acgt_counts[idx] += 1
            node = name_to_node[name]
            nodes.extend([node]*4)
            sequences.append([nt_to_state[c] for c in sequence])
            variables.extend([0, 1, 2, 3])
    # arr[node, site, variable]
    arr = np.array(sequences, dtype=int)
    arr = np.transpose(arr, (1, 0, 2))
    iid_observations = np.reshape(arr, (arr.shape[0], -1))
    iid_observations = iid_observations.tolist()
    pi = (acgt_counts / acgt_counts.sum()).tolist()
    return nodes, variables, iid_observations, pi


def get_process_definition(pi, kappa):
    # Get the structure of the rate matrix.
    ts_info, tv_info = get_transition_components(pi)
    (
            ts_row_states,
            ts_column_states,
            ts_row_probabilities,
            ts_coefficients) = ts_info
    (
            tv_row_states,
            tv_column_states,
            tv_row_probabilities,
            tv_coefficients) = tv_info

    # Put together the rate matrix.
    row_states = ts_row_states + tv_row_states
    column_states = ts_column_states + tv_column_states
    row_probabilities = ts_row_probabilities + tv_row_probabilities
    coefficients = ts_coefficients + tv_coefficients
    nts = len(ts_row_states)
    ntv = len(tv_row_states)
    rates = [c * r for c, r in zip(coefficients, [kappa]*nts + [1]*ntv)]
    expected_rate = np.dot(row_probabilities, rates)
    rates = [r / expected_rate for r in rates]

    # Assemble and return the process definition.
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = rates)
    return process_definition


def get_transition_components(pi):
    # Returns ts_info, tv_info,
    # each of which has (row_states, column_states, coefficients).
    ts_row_states = []
    ts_column_states = []
    ts_row_probabilities = []
    ts_coefficients = []
    tv_row_probabilities = []
    tv_row_states = []
    tv_column_states = []
    tv_coefficients = []
    ident = np.identity(4, dtype=int).tolist()
    # ts: A<->G, C<->T
    ts = ((0, 2), (2, 0), (1, 3), (3, 1))
    for i in range(4):
        for j in range(4):
            if i != j:
                rs = ident[i]
                cs = ident[j]
                if (i, j) in ts:
                    ts_row_states.append(rs)
                    ts_column_states.append(cs)
                    ts_row_probabilities.append(pi[i])
                    ts_coefficients.append(pi[j])
                else:
                    tv_row_states.append(rs)
                    tv_column_states.append(cs)
                    tv_row_probabilities.append(pi[i])
                    tv_coefficients.append(pi[j])
    ts_info = (
            ts_row_states,
            ts_column_states,
            ts_row_probabilities,
            ts_coefficients)
    tv_info = (
            tv_row_states,
            tv_column_states,
            tv_row_probabilities,
            tv_coefficients)
    return ts_info, tv_info


def main():

    # Read the tree.
    name_to_node, edges = get_tree_info()
    edge_count = len(edges)
    node_count = edge_count + 1

    # Read the alignment.
    info = get_nucleotide_alignment_info(name_to_node)
    nodes, variables, iid_observations, pi = info
    nsites = len(iid_observations)

    # Define the tree component of the scene
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = [0.01] * edge_count,
            edge_processes = [0] * edge_count)

    # Get the structure of the rate matrix.
    ts_info, tv_info = get_transition_components(pi)
    (
            ts_row_states,
            ts_column_states,
            ts_row_probabilities,
            ts_coefficients) = ts_info
    (
            tv_row_states,
            tv_column_states,
            tv_row_probabilities,
            tv_coefficients) = tv_info

    row_states = ts_row_states + tv_row_states
    column_states = ts_column_states + tv_column_states

    # Define the root distribution.
    ident = np.identity(4, dtype=int).tolist()
    root_prior = dict(
            states = ident,
            probabilities = [0.25, 0.25, 0.25, 0.25])

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

    # Define some requests.
    # These include the log likelihood
    # and the sum of transition count expectations.
    requests = [
            {"property" : "SNNLOGL"},
            {
                "property" : "SDWDWEL",
                "state_reduction" : {
                    "states" : np.identity(4).tolist(),
                    "weights" : [1] * len(row_states)
                }
            },
            {
                "property" : "SSNTRAN",
                "transition_reduction" : {
                    "row_states" : ts_row_states,
                    "column_states" : ts_column_states,
                    "weights" : [1] * len(ts_row_states)
                }
            },
            {
                "property" : "SSNTRAN",
                "transition_reduction" : {
                    "row_states" : tv_row_states,
                    "column_states" : tv_column_states,
                    "weights" : [1] * len(tv_row_states)
                }
            }]

    # Request some stuff.
    j_in = {
            "scene" : scene,
            "requests" : requests
            }

    arr = []
    j_out = None
    kappa = 2.0
    for i in range(8):

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
            j_in['scene']['tree']['edge_rate_scaling_factors'] = edge_rates

        process_definition = get_process_definition(pi, kappa)
        j_in['scene']['process_definitions'] = [process_definition]
        j_out = jsonctmctree.interface.process_json_in(j_in)
        arr.append(copy.deepcopy(j_out))

    print(json.dumps(arr, indent=4))


main()
