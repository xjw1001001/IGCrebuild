from __future__ import print_function, division, absolute_import

import functools
import json

import numpy as np
from numpy.testing import assert_equal
from scipy.misc import logsumexp
from scipy.optimize import minimize
import pyparsing

import jsonctmctree.interface
from jsonctmctree.extras import optimize_em

from modelutil import pack, unpack


s_tree = """((((((cerevisiae,paradoxus),mikatae),kudriavzevii),
bayanus),castellii),kluyveri)"""

FASTA_FILENAME = 'YDR502C_YLR180W_input.fasta'
#FASTA_FILENAME = 'yeast.paralogs.fasta'

PARALOGS = ('YDR502C', 'YLR180W')
#PARALOGS = ('YAL056W', 'YOR371C')


def custom_pack(distn, kappa, tau, rates):
    return pack(distn, kappa, tau, hardcoded_rate_deflate(rates))

def custom_unpack(rate_expansion, X):
    distn, kappa, tau, rates = unpack(X)
    return distn, kappa, tau, hardcoded_rate_expand(rate_expansion, rates)

def hardcoded_rate_deflate(rates):
    new_rates = np.concatenate((rates[:1], rates[2:]))
    assert_equal(len(new_rates) + 1, len(rates))
    return new_rates

def hardcoded_rate_expand(outgroup_rate, rates):
    new_rates = np.concatenate((rates[:1], [outgroup_rate], rates[1:]))
    assert_equal(len(new_rates), len(rates) + 1)
    return new_rates


def get_hardcoded_tree_and_rates(outgroup_rate):
    # This is for compatibility with Xiang's implementation.
    # It is an alternative to build_tree.
    names = [
            'N0', 'N1', 'N2', 'N3', 'N4', 'N5',
            'kluyveri', 'castellii', 'bayanus', 'kudriavzevii',
            'mikatae', 'cerevisiae', 'paradoxus']
    name_to_node = {x : i for i, x in enumerate(names)}
    name_edges = [
            ('N0','N1'),
            ('N0','kluyveri'),
            ('N1','N2'),
            ('N1','castellii'),
            ('N2','N3'),
            ('N2','bayanus'),
            ('N3','N4'),
            ('N3','kudriavzevii'),
            ('N4','N5'),
            ('N4','mikatae'),
            ('N5','cerevisiae'),
            ('N5','paradoxus'),
            ]
    edges = [[name_to_node[a], name_to_node[b]] for a, b in name_edges]
    rates = [
            1.323666e-01,

            1.141184e-05,
            #outgroup_rate,

            3.208034e-02,
            1.167843e-01,
            2.156341e-02,
            2.942265e-02,
            1.778006e-02,
            4.161289e-02,
            1.778255e-02,
            4.193464e-02,
            2.272669e-02,
            2.114252e-02,
            ]
    return name_to_node, edges, rates


def build_tree(parent, root, node, name_to_node, edges):
    if parent is not None:
        edges.append((parent, node))
    neo = node + 1
    if isinstance(root, basestring):
        name_to_node[root] = node
    else:
        for element in root:
            neo = build_tree(node, element, neo, name_to_node, edges)
    return neo

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

def objective_and_gradient(scene, outgroup_length, X):
    #print('parameter estimates passed to objective and gradient calculator:')
    #print(X)
    #print(np.exp(X))
    print_estimates(outgroup_length, X)

    delta = 1e-8
    distn, kappa, tau, rates = custom_unpack(outgroup_length, X)
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
    j_out = jsonctmctree.interface.process_json_in(j_in,
            #debug=True,
            )
    status = j_out['status']
    assert_equal(status, 'feasible')
    log_likelihood, edge_gradient = j_out['responses']
    cost = -log_likelihood

    # For each non-edge-specific parameter get finite-differences
    # approximation of the gradient.
    nedges = len(scene['tree']['row_nodes'])
    nparams = len(X) - nedges
    gradient = []
    for i in range(nparams):
        W = np.copy(X)
        W[i] += delta
        distn, kappa, tau, rates = custom_unpack(outgroup_length, W)
        process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
        scene['root_prior'] = root_prior
        scene['process_definitions'] = [process_defn]
        j_in = {
                'scene' : scene,
                'requests' : [log_likelihood_request]
                }
        j_out = jsonctmctree.interface.process_json_in(j_in)
        ll = j_out['responses'][0]
        c = -ll
        slope = (c - cost) / delta
        gradient.append(slope)
    gradient.extend([-x for x in edge_gradient])
    gradient = np.array(gradient)

    # Return cost and gradient.
    print(cost)
    print('gradient:')
    print(gradient)
    print()
    return cost, gradient


def compute_log_likelihood_for_one_branch_length(outgroup_length):
    # Flatten the tree into a list of node indices and a list of edges.
    #tree = s_tree.replace(',', ' ')
    #nestedItems = pyparsing.nestedExpr(opener='(', closer=')')
    #tree = (nestedItems + pyparsing.stringEnd).parseString(tree).asList()[0]
    #name_to_node = {}
    #edges = []
    #build_tree(None, tree, 0, name_to_node, edges)
    name_to_node, edges, rates = get_hardcoded_tree_and_rates(outgroup_length)

    # Spam the name and node and edges.
    for name, node in name_to_node.items():
        print(node, ':', name)
    for edge in edges:
        print(edge)
    print()

    # Define suffixes indicating paralogs.
    paralog_to_variable = {
            PARALOGS[0].lower() : 0,
            PARALOGS[1].lower() : 1}

    # Read the data (in this case, alignment columns).
    nodes = []
    variables = []
    rows = []
    with open(FASTA_FILENAME) as fin:
        while True:
            line = fin.readline().strip().lower()
            if not line:
                break
            name = line[1:-7]
            paralog = line[-7:]
            seq = fin.readline().strip()
            row = ['ACGT'.index(x) for x in seq]
            nodes.append(name_to_node[name])
            variables.append(paralog_to_variable[paralog])
            rows.append(row)
    columns = [list(x) for x in zip(*rows)]

    # Interpret the tree.
    observed_node_count = len(nodes)
    edge_count = len(edges)
    row_nodes, column_nodes = zip(*edges)

    print('number of sites in the alignment:', len(columns))
    print('number of sequences:', observed_node_count)

    # Compute the empirical distribution of the nucleotides.
    counts = np.zeros(4)
    for k in np.ravel(columns):
        counts[k] += 1
    empirical_pi = counts / counts.sum()

    distn = empirical_pi
    #distn = np.ones(4) / 4
    #distn = np.array([2.5125e-01, 2.3797e-01, 2.0262e-01, 3.0815e-01])

    #kappa = 5.944071e+00
    #tau = 3.140967e+00

    # This has been replaced by hardcoded rates.
    #rates = [0.1] * edge_count

    kappa = 6.0
    tau = 3.0
    process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
    scene = {
            "node_count" : edge_count + 1,
            "process_count" : 1,
            "state_space_shape" : [4, 4],
            "tree" : {
                "row_nodes" : list(row_nodes),
                "column_nodes" : list(column_nodes),
                "edge_rate_scaling_factors" : rates,
                "edge_processes" : [0] * edge_count
                },
            "root_prior" : root_prior,
            "process_definitions" : [process_defn],
            "observed_data" : {
                "nodes" : nodes,
                "variables" : variables,
                "iid_observations" : columns
                }
            }

    # Update the edge rates using 5 iterations of EM.
    #rates = optimize_em(scene, None, 5)
    #assert_equal(len(rates), len(edges))

    X = custom_pack(distn, kappa, tau, rates)
    f = functools.partial(objective_and_gradient, scene, outgroup_length)
    result = minimize(f, X, jac=True, method='L-BFGS-B')
    print('post-search value of objective function:', result.fun)
    X = result.x
    print('post-search results:')
    print_estimates(outgroup_length, X)
    print()

    # Update the edge rates using 5 iterations of EM.
    # Then re-run the optimization.
    """
    distn, kappa, tau, rates = custom_unpack(outgroup_length, X)
    process_defn, root_prior = get_process_defn_and_prior(distn, kappa, tau)
    scene['root_prior'] = root_prior
    scene['process_definitions'] = [process_defn]
    rates = optimize_em(scene, None, 5)
    X = custom_pack(distn, kappa, tau, rates)
    f = functools.partial(objective_and_gradient, scene, outgroup_length)
    result = minimize(f, X, jac=True, method='L-BFGS-B')
    print('final results:')
    print('objective function value:', result.fun)
    X = result.x
    print_estimates(outgroup_length, X)
    print()
    """

    return result.fun


def print_estimates(outgroup_length, X):
    distn, kappa, tau, rates = custom_unpack(outgroup_length, X)
    print('nucleotide distribution:')
    for nt, p in zip('ACGT', distn):
        print('  ', nt, ':', p)
    print('kappa:', kappa)
    print('tau:', tau)
    print('edge rate scaling factors:')
    for r in rates:
        print('  ', r)


def main():
    np.set_printoptions(threshold=100000)

    log_likelihoods = []
    nsamples = 15
    log_outgroup_rates = np.linspace(-10, -3, num=nsamples)
    outgroup_rates = np.exp(log_outgroup_rates)
    for outgroup_rate in outgroup_rates:
        ll = compute_log_likelihood_for_one_branch_length(outgroup_rate)
        log_likelihoods.append(ll)
        print(ll)

    print('outgroup rate and maximum log likelihood:')
    for outgroup_rate, ll in zip(outgroup_rates, log_likelihoods):
        print(outgroup_rate, ll, sep='\t')

    print('draw the plot:')
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.plot(outgroup_rates, log_likelihoods, marker='o', linestyle='--')
    plt.xscale('log')
    plt.xlabel('outgroup branch length')
    plt.ylabel('neg log likelihood')
    plt.title('outgroup branch discretization for HKY+IGC YDR502C YLR180W')
    plt.savefig('discretized-outgroup-plot.svg')


main()
