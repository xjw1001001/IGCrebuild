"""
Use the 'extras' module to help with the inference.

This script has no clock-like constraint on branch length,
but it has a zero inter-locus gene-conversion constraint.
It uses EM only for an improvement of the initial guesses of the
edge-specific rate scaling factors.

The 'extras' module in the jsonctmctree Python package
is used for some EM and quasi-newton optimization.

"""
from __future__ import print_function, division

import argparse
import itertools

import pyparsing

import numpy as np
from numpy.testing import assert_
import scipy.special

from jsonctmctree.interface import process_json_in
from jsonctmctree.extras import optimize_em, optimize_quasi_newton


###############################################################################
# Stuff in this section is related to the data.


def gen_paragraphs(fin):
    lines = []
    for line in fin:
        line = line.rstrip()
        if line:
            lines.append(line)
        else:
            if lines:
                yield lines
                lines = []
    if lines:
        yield lines


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


def get_tree_info(tree_string):
    # Return a dictionary mapping name to node index,
    # and return a list of edges as ordered pairs of node indices.
    assert_(tree_string.endswith(';'))
    tree = tree_string[:-1].replace(',', ' ')
    nestedItems = pyparsing.nestedExpr(opener='(', closer=')')
    tree = (nestedItems + pyparsing.stringEnd).parseString(tree).asList()[0]
    name_to_node = {}
    edges = []
    _help_build_tree(None, tree, 0, name_to_node, edges)
    return name_to_node, edges


def parse_full_name(full_name, paralog_names):
    # Return (species_name, paralog_name_index).
    for i, paralog_name in enumerate(paralog_names):
        if full_name.endswith(paralog_name):
            species_name = full_name[:-len(paralog_name)]
            return species_name, i
    raise Exception(full_name)


def get_alignment_info(fasta_fd, name_to_node, paralog_names):
    """
    Read the alignment data.

    Parameters
    ----------
    fasta_fd : open file-like object
        The nucleotide alignment.
    name_to_node : dict
        Map the species name to the tree node.
    paralog_names : sequence of strings
        Sequence of paralog names.

    Returns
    -------
    nodes : sequence of integers
        Sequence of observable nodes.
        Nodes may be repeated if multiple variables are observable per node.
    variables : sequence of integers
        Sequence of observable variables.
    columns : sequence of integer lists
        Sequence of observation lists.

    """
    nodes = []
    variables = []
    rows = []
    for lines in gen_paragraphs(fasta_fd):
        if len(lines) != 2:
            raise Exception('expected two lines per paragraph')

        # Process the name line, containing the species and paralog.
        name_line = lines[0].strip()
        name, variable = parse_full_name(name_line, paralog_names)
        nodes.append(name_to_node[name])
        variables.append(variable)

        # Process the sequence line, containing the DNA sequence.
        sequence_line = lines[1].strip()
        row = ['ACGT'.index(x) for x in sequence_line]
        rows.append(row)

    columns = [list(x) for x in zip(*rows)]
    return nodes, variables, columns


###############################################################################
# Stuff in this section is specific to the parameterization details
# of the model.


def pack_global_params(pi, kappa):
    a, c, g, t = pi
    acgt = pi.sum()
    at = a+t
    cg = c+g
    a_div_at = a / at
    c_div_cg = c / cg
    arr = np.concatenate([
        scipy.special.logit([at, a_div_at, c_div_cg]),
        np.log([kappa])])
    return arr


def unpack_global_params(P):
    nt_info = scipy.special.expit(P[0:3])
    at, a_div_at, c_div_cg = nt_info
    a = a_div_at * at
    t = (1 - a_div_at) * at
    cg = 1 - at
    c = c_div_cg * cg
    g = (1 - c_div_cg) * cg
    pi = np.array([a, c, g, t])
    kappa = np.exp(P[-1])
    return pi, kappa


def _get_process_definitions(P):
    # This is called within the optimization.
    pi, kappa = unpack_global_params(P)
    return [get_joint_hky_process_definition(pi, kappa)]


def _get_root_prior(P):
    # This is called within the optimization.
    pi, kappa = unpack_global_params(P)
    return get_root_prior(pi)


###############################################################################
# Stuff in this section is specific to the joint paralog HKY model
# and does not care about the input formats of the data.


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


def gen_joint_hky():
    # Yield a tuple for each state transition.
    # The tuple consists of the multivariate initial state,
    # the multivariate final state,
    # a ts indicator, a tv indicator,
    # and a final mutational nucleotide index.

    # Precompute the transitions out of each nucleotide state.
    row_idx_to_info = [[] for i in range(4)]
    for info in gen_hky():
        i, j, ts, tv = info
        row_idx_to_info[i].append(info)

    # Compute the joint state transitions.
    for ia, ib in itertools.product(range(4), repeat=2):

        # Iterate over all transitions for the first nucleotide.
        for i, j, ts, tv in row_idx_to_info[ia]:
            ja, jb = j, ib
            yield [ia, ib], [ja, jb], ts, tv, j

        # Iterate over all transitions for the second nucleotide.
        for i, j, ts, tv in row_idx_to_info[ib]:
            ja, jb = ia, j
            yield [ia, ib], [ja, jb], ts, tv, j


def get_joint_hky_process_definition(pi, kappa):
    # Note that the expected rate normalization is for
    # only a single site, not for both sites.
    # This is intentional.
    # So the expected number of changes along an edge
    # is about twice the edge rate scaling factor of that edge.
    expected_rate = get_expected_univariate_rate(pi, kappa)
    info = get_unnormalized_transitions(pi, kappa)
    row_states, column_states, transition_rates = info
    normalized_rates = [r / expected_rate for r in transition_rates]
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = normalized_rates)
    return process_definition


def get_root_prior(pi):
    root_prior = dict(
            states = [[i, i] for i in range(4)],
            probabilities = list(pi))
    return root_prior


def get_unnormalized_ts_transitions(pi, kappa):
    row_states = []
    column_states = []
    transition_rates = []
    for row_state, column_state, ts, tv, j in gen_joint_hky():
        if ts:
            exit_rate = (kappa * ts + tv) * pi[j]
            row_states.append(row_state)
            column_states.append(column_state)
            transition_rates.append(exit_rate)
    return row_states, column_states, transition_rates


def get_unnormalized_tv_transitions(pi, kappa):
    row_states = []
    column_states = []
    transition_rates = []
    for row_state, column_state, ts, tv, j in gen_joint_hky():
        if tv:
            exit_rate = (kappa * ts + tv) * pi[j]
            row_states.append(row_state)
            column_states.append(column_state)
            transition_rates.append(exit_rate)
    return row_states, column_states, transition_rates


def get_unnormalized_transitions(pi, kappa):
    ts_row, ts_col, ts_rate = get_unnormalized_ts_transitions(pi, kappa)
    tv_row, tv_col, tv_rate = get_unnormalized_tv_transitions(pi, kappa)
    row_states = ts_row + tv_row
    column_states = ts_col + tv_col
    transition_rates = ts_rate + tv_rate
    return row_states, column_states, transition_rates


def get_expected_univariate_rate(pi, kappa):
    raw_exit_rates = get_unnormalized_univariate_exit_rates(pi, kappa)
    return np.dot(pi, raw_exit_rates)


def get_unnormalized_univariate_exit_rates(pi, kappa):
    exit_rates = [0, 0, 0, 0]
    for i, j, ts, tv in gen_hky():
        rate = (kappa * ts + tv) * pi[j]
        exit_rates[i] += rate
    return exit_rates


###############################################################################
# The main script function.


def main(args):

    # Get the paralog names.
    paralog_names = args.paralogs

    # Read the tree.
    with open(args.tree) as fin:
        tree_string = fin.read().strip()
    name_to_node, edges = get_tree_info(tree_string)
    edge_count = len(edges)
    node_count = edge_count + 1

    # Read the alignment.
    with open(args.alignment) as alignment_fd:
        info = get_alignment_info(alignment_fd, name_to_node, paralog_names)
    nodes, variables, iid_observations = info
    nsites = len(iid_observations)

    print('number of sites in the alignment:', nsites)
    print('number of sequences:', len(nodes))

    # Compute the empirical distribution of the nucleotides.
    counts = np.zeros(4)
    for k in np.ravel(iid_observations):
        counts[k] += 1
    empirical_pi = counts / counts.sum()

    # Initialize some guesses.
    edge_rates = [0.01] * edge_count
    pi = empirical_pi
    kappa = 2.0

    # Define the tree component of the scene
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0] * edge_count)

    # Define the root distribution.
    root_prior = get_root_prior(pi)

    # Define the observed data.
    observed_data = dict(
            nodes = nodes,
            variables = variables,
            iid_observations = iid_observations)

    # Assemble the scene.
    process_defn = get_joint_hky_process_definition(pi, kappa)
    scene = dict(
            node_count = node_count,
            process_count = 1,
            state_space_shape = [4, 4],
            tree = tree,
            root_prior = root_prior,
            process_definitions = [process_defn],
            observed_data = observed_data)

    print('computing the log likelihood...')

    # Ask for the log likelihood, summed over sites.
    log_likelihood_request = dict(property = 'SNNLOGL')
    j_in = dict(
            scene = scene,
            requests = [log_likelihood_request])
    j_out = process_json_in(j_in)
    print(j_out)

    print('updating edge specific rate scaling factors using EM...')

    # Use the generic EM edge rate scaling factor updating function.
    observation_reduction = None
    em_iterations = 1
    edge_rates = optimize_em(scene, observation_reduction, em_iterations)

    # Update the scene to reflect the edge rates.
    print('updated edge rate scaling factors:')
    print(edge_rates)
    scene['tree']['edge_rate_scaling_factors'] = edge_rates

    print('checking log likelihood after having updated edge rates...')

    # Check the log likelihood again.
    j_in = dict(
            scene = scene,
            requests = [log_likelihood_request])
    j_out = process_json_in(j_in)
    print(j_out)

    print('computing the maximum likelihood estimates...')

    # Improve the estimates using a numerical search.
    P0 = pack_global_params(pi, kappa)
    B0 = np.log(edge_rates)
    verbose = False
    observation_reduction = None
    result, P_opt, B_opt = optimize_quasi_newton(
            verbose,
            scene,
            observation_reduction,
            _get_process_definitions,
            _get_root_prior,
            P0, B0)

    # Unpack and report the results.
    pi, kappa = unpack_global_params(P_opt)
    edge_rates = np.exp(B_opt)
    print('negative log likelihood:', result.fun)
    print('nucleotide distribution:')
    for nt, p in zip('ACGT', pi):
        print(nt, ':', p)
    print('kappa:', kappa)
    print('edge rates:')
    print(edge_rates)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment', required=True,
            help='alignment file')
    parser.add_argument('--tree', required=True,
            help='tree file')
    parser.add_argument('--paralogs', nargs='+', required=True,
            help='paralog names')
    args = parser.parse_args()
    main(args)
