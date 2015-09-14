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
    sequences = []
    nodes = []
    variables = []
    with open('vegetables.rbcL.txt') as fin:
        lines = fin.readlines()
        header = lines[0]
        for line in lines[1:]:
            name, sequence = line.strip().split()
            node = name_to_node[name]
            nodes.extend([node]*4)
            sequences.append([nt_to_state[c] for c in sequence])
            variables.extend([0, 1, 2, 3])
    # arr[node, site, variable]
    arr = np.array(sequences, dtype=int)
    arr = np.transpose(arr, (1, 0, 2))
    iid_observations = np.reshape(arr, (arr.shape[0], -1))
    iid_observations = iid_observations.tolist()
    return nodes, variables, iid_observations


def get_process_definition(kappa):
    # Get the structure of the rate matrix.
    info = get_transition_components()
    ts_row_states, ts_column_states, tv_row_states, tv_column_states = info

    # Put together the rate matrix.
    rate = 2 + kappa
    row_states = ts_row_states + tv_row_states
    column_states = ts_column_states + tv_column_states
    rates = [kappa/rate]*len(ts_row_states) + [1/rate]*len(tv_row_states)

    # Assemble and return the process definition.
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = rates)
    return process_definition


def get_transition_components():
    # Return ts_row_states, ts_col_states, tv_row_states, tv_col_states
    ts_row_states = []
    ts_column_states = []
    tv_row_states = []
    tv_column_states = []
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
                else:
                    tv_row_states.append(rs)
                    tv_column_states.append(cs)
    return ts_row_states, ts_column_states, tv_row_states, tv_column_states


def main():

    # Read the tree.
    name_to_node, edges = get_tree_info()
    edge_count = len(edges)
    node_count = edge_count + 1

    # Read the alignment.
    info = get_nucleotide_alignment_info(name_to_node)
    nodes, variables, iid_observations = info
    nsites = len(iid_observations)

    # Define the tree component of the scene
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = [0.01] * edge_count,
            edge_processes = [0] * edge_count)

    # Get the structure of the rate matrix.
    info = get_transition_components()
    ts_row_states, ts_column_states, tv_row_states, tv_column_states = info
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
                "property" : "SDNTRAN",
                "transition_reduction" : {
                    "row_states" : row_states,
                    "column_states" : column_states,
                    "weights" : [1 / nsites] * len(row_states)
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
            ll, edge_rates, ts, tv = j_out['responses']
            kappa = 2 * (ts / tv)
            j_in['scene']['tree']['edge_rate_scaling_factors'] = edge_rates

        j_in['scene']['process_definitions'] = [get_process_definition(kappa)]
        j_out = jsonctmctree.interface.process_json_in(j_in)
        arr.append(copy.deepcopy(j_out))

    print(json.dumps(arr, indent=4))
    print('kappa estimate:', kappa)


main()
