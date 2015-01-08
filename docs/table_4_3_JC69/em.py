from __future__ import print_function, division, absolute_import

import copy
import json

import numpy as np
from numpy.testing import assert_equal

from scipy.misc import logsumexp
from scipy.optimize import minimize
import pyparsing

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
    print(arr.shape)
    arr = np.transpose(arr, (1, 0, 2))
    iid_observations = np.reshape(arr, (arr.shape[0], -1))
    print(iid_observations.shape)
    iid_observations = iid_observations.tolist()
    #iid_observations = [list(x) for x in zip(*sequences)]
    #iid_observations = [x for x in iid_observations if -1 not in x]
    print(nodes)
    print(variables)
    #print(iid_observations[0])
    return nodes, variables, iid_observations


def main():

    # Read the tree.
    name_to_node, edges = get_tree_info()
    edge_count = len(edges)
    node_count = edge_count + 1

    # Read the alignment.
    info = get_nucleotide_alignment_info(name_to_node)
    nodes, variables, iid_observations = info

    # Define the tree component of the scene
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = [0.01] * len(edges),
            edge_processes = [0] * len(edges))

    # Define the Jukes-Cantor process.
    # Rates are scaled so that the exit rate is 1 from every state.
    row_states = []
    column_states = []
    rates = []
    for i in range(4):
        for j in range(4):
            if i != j:
                u = [0]*4
                u[i] = 1
                v = [0]*4
                v[j] = 1
                row_states.append(u)
                column_states.append(v)
                rates.append(1/3)
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = rates)

    # Define the root distribution.
    root_prior = dict(
            states = [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]],
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
            process_definitions = [process_definition],
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
                    "weights" : [1] * len(row_states)
                }
            }]

    # Request some stuff.
    j_in = {
            "scene" : scene,
            "requests" : requests
            }

    arr = []
    j_out = None
    nsites = len(iid_observations)
    print('nsites:', nsites)
    print('edge count:', edge_count)
    for i in range(10):
        if j_out is None:
            j_out = jsonctmctree.interface.process_json_in(j_in)
        else:
            ll, trans_counts = j_out['responses']
            edge_rates = [t/nsites for t in trans_counts]
            j_in['scene']['tree']['edge_rate_scaling_factors'] = edge_rates
            j_out = jsonctmctree.interface.process_json_in(j_in)
            print(ll)
            print(edge_rates)
        arr.append(copy.deepcopy(j_out))
        #print(j_out)

    #print(json.dumps(arr, indent=4))


main()
