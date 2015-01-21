from __future__ import print_function, division

from functools import partial
import itertools
from StringIO import StringIO
import copy

import numpy as np
from numpy.testing import assert_equal, assert_
import scipy.optimize

import dendropy

import jsonctmctree.interface


_code = """
0	ala	gct
1	ala	gcc
2	ala	gca
3	ala	gcg
4	arg	cgt
5	arg	cgc
6	arg	cga
7	arg	cgg
8	arg	aga
9	arg	agg
10	asn	aat
11	asn	aac
12	asp	gat
13	asp	gac
14	cys	tgt
15	cys	tgc
16	gln	caa
17	gln	cag
18	glu	gaa
19	glu	gag
20	gly	ggt
21	gly	ggc
22	gly	gga
23	gly	ggg
24	his	cat
25	his	cac
26	ile	att
27	ile	atc
28	ile	ata
29	leu	tta
30	leu	ttg
31	leu	ctt
32	leu	ctc
33	leu	cta
34	leu	ctg
35	lys	aaa
36	lys	aag
37	met	atg
38	phe	ttt
39	phe	ttc
40	pro	cct
41	pro	ccc
42	pro	cca
43	pro	ccg
44	ser	tct
45	ser	tcc
46	ser	tca
47	ser	tcg
48	ser	agt
49	ser	agc
50	thr	act
51	thr	acc
52	thr	aca
53	thr	acg
54	trp	tgg
55	tyr	tat
56	tyr	tac
57	val	gtt
58	val	gtc
59	val	gta
60	val	gtg
61	stop	taa
62	stop	tga
63	stop	tag
"""

# This is an official python itertools recipe.
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.izip_longest(*args, fillvalue=fillvalue)


def get_deltas(a, b):
    deltas = []
    for x, y in zip(a, b):
        if x != y:
            deltas.append((x, y))
    return deltas


def pack_acgt(pi):
    a, c, g, t = pi
    ag = a+g  # purines
    ct = c+t  # pyrimidines
    a_div_ag = a / ag
    c_div_ct = c / ct
    return scipy.special.logit([ag, a_div_ag, c_div_ct])


def unpack_acgt(packed_acgt):
    ag, a_div_ag, c_div_ct = scipy.special.expit(packed_acgt)
    ct = 1 - ag
    a = a_div_ag * ag
    g = ag - a
    c = c_div_ct * ct
    t = ct - c
    return np.array([a, c, g, t])


def read_newick(fin):
    # use dendropy to read this newick file
    t = dendropy.Tree(stream=fin, schema='newick')
    nodes = list(t.preorder_node_iter())
    id_to_idx = {id(n) : i for i, n in enumerate(nodes)}
    edges = []
    edge_rates = []
    for dendro_edge in t.preorder_edge_iter():
        if dendro_edge.tail_node and dendro_edge.head_node:
            na = id_to_idx[id(dendro_edge.tail_node)]
            nb = id_to_idx[id(dendro_edge.head_node)]
            edges.append((na, nb))
            edge_rates.append(dendro_edge.length)
    name_to_node = {str(n.taxon) : id_to_idx[id(n)] for n in t.leaf_nodes()}
    return edges, edge_rates, name_to_node


def codon_to_triple(codon):
    return ['ACGT'.index(c) for c in codon]


def gen_mg94_structure_by_rows(codon_residue_pairs):
    transitions = (('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C'))
    for i, (ca, ra) in enumerate(codon_residue_pairs):
        info_tuples = []
        for j, (cb, rb) in enumerate(codon_residue_pairs):
            if i != j:
                deltas = get_deltas(ca, cb)
                if len(deltas) == 1:
                    delta = deltas[0]
                    if delta in transitions:
                        ts, tv = 1, 0
                    else:
                        ts, tv = 0, 1
                    if ra == rb:
                        non, syn = 0, 1
                    else:
                        non, syn = 1, 0
                    nt = 'ACGT'.index(delta[1])
                    info = (i, j, ts, tv, non, syn, nt)
                    info_tuples.append(info)
        if info_tuples:
            yield i, info_tuples


def gen_mg94_structure(codon_residue_pairs):
    # Yield (i, j, ts, tv, non, syn, nt).
    for ia, infos in gen_mg94_structure_by_rows(codon_residue_pairs):
        for info in infos:
            yield info


def gen_mg94_exit_rates(pi, kappa, omega, codon_residue_pairs):
    # This does not include normalization.
    for ia, infos in gen_mg94_structure_by_rows(codon_residue_pairs):
        exit_rate = 0
        for i, j, ts, tv, non, syn, nt in infos:
            assert_equal(ia, i)
            rate = (kappa * ts + tv) * (omega * non + syn) * pi[nt]
            exit_rate += rate
        yield exit_rate


def gen_geneconv_structure_by_rows(pi, kappa, omega, tau, codon_residue_pairs):

    # Define the MG94 structure by rows.
    # This will be used as an ingredient
    # in building the geneconv rate matrix structure.
    mg94_by_rows = dict(gen_mg94_structure_by_rows(codon_residue_pairs))

    # Iterate over initial codon pairs.
    # These are the multivariate 'row states';
    # for each such multivariate row state,
    # we will yield the multivariate row state followed by a sequence
    # of transition information tuples.
    for ia, ib in itertools.product(range(61), repeat=2):
        multivariate_row_state = [ia, ib]
        info_tuples = []

        # Iterate over all substitution transitions for the first codon.
        # Include the interlocus gene conversion if appropriate.
        homogenized_to_site_b = False
        for i, j, ts, tv, non, syn, nt in mg94_by_rows[ia]:
            assert_equal(ia, i)
            ja, jb = j, ib
            if ja == jb:
                assert_(not homogenized_to_site_b)
                homogenized_to_site_b = True
                hom = 1
            else:
                hom = 0
            multi_info = [ia, ib], [ja, jb], ts, tv, non, syn, nt, hom
            info_tuples.append(multi_info)

        # Iterate over all substitution transitions for the second codon.
        # Include the interlocus gene conversion if appropriate.
        homogenized_to_site_a = False
        for i, j, ts, tv, non, syn, nt in mg94_by_rows[ib]:
            assert_equal(ib, i)
            ja, jb = ia, j
            if ja == jb:
                assert_(not homogenized_to_site_a)
                homogenized_to_site_a = True
                hom = 1
            else:
                hom = 0
            multi_info = [ia, ib], [ja, jb], ts, tv, non, syn, nt, hom
            info_tuples.append(multi_info)

        # If only one of the homogenizations has occurred,
        # then this would indicate a bug.
        if homogenized_to_site_a and not homogenized_to_site_b:
            raise Exception
        if homogenized_to_site_b and not homogenized_to_site_a:
            raise Exception

        # If the initial pair of codons are different from each other,
        # and if homogenization has not yet occurred,
        # then add one more interlocus gene conversion in each direction.
        if ia != ib:
            if not homogenized_to_site_a and not homogenized_to_site_b:

                # The interlocus gene conversion rate will
                # depend on whether the residues are the same or not.
                ia_codon, ia_residue = codon_residue_pairs[ia]
                ib_codon, ib_residue = codon_residue_pairs[ib]

                # This is neither a transition nor a transversion.
                ts = 0
                tv = 0
                syn = 1 if ia_residue == ib_residue else 0
                non = 1 - syn
                hom = 1
                nt = None

                # Add interlocus gene conversion in one direction.
                multi_info = [ia, ib], [ib, ib], ts, tv, non, syn, nt, hom
                info_tuples.append(multi_info)

                # Add interlocus gene conversion in the other direction.
                multi_info = [ia, ib], [ia, ia], ts, tv, non, syn, nt, hom
                info_tuples.append(multi_info)

        # This is a prototype.
        yield multivariate_row_state, info_tuples


def get_geneconv_process_definition(
        pi, kappa, omega, tau, codon_distribution, codon_residue_pairs):

    # First, compute an expected rate for the univariate MG94 process.
    mg94_exit_rates = list(gen_mg94_exit_rates(
        pi, kappa, omega, codon_residue_pairs))
    mg94_expected_rate = np.dot(codon_distribution, mg94_exit_rates)

    # Next iterate over the inter-locus gene conversion process structure.
    # For each transition, compute the rate.
    row_states = []
    column_states = []
    transition_rates = []
    for initial_state, infos in gen_geneconv_structure_by_rows(
            pi, kappa, omega, tau, codon_residue_pairs):
        for row_state, column_state, ts, tv, non, syn, nt, hom in infos:
            ia, ib = row_state
            ja, jb = column_state
            assert_equal(initial_state, row_state)
            rmut = 0
            rhom = 0
            if nt is not None:
                rmut = (kappa * ts + tv) * (omega * non + syn) * pi[nt]
            if hom:
                rhom = (omega * non + syn) * tau
            transition_rate = (rmut + rhom) / mg94_expected_rate
            row_states.append(row_state)
            column_states.append(column_state)
            transition_rates.append(transition_rate)

    # Assemble the process definition.
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates)
    return process_definition


def pack(pi, kappa, omega, tau, edge_rates):
    return np.concatenate([
        pack_acgt(pi),
        np.log([kappa, omega, tau]),
        np.log(edge_rates)])


def unpack(X):
    pi = unpack_acgt(X[0:3])
    misc_params = np.exp(X[3:3+3])
    kappa, omega, tau = misc_params
    edge_rates = np.exp(X[6:])
    return pi, kappa, omega, tau, edge_rates


def print_packed(X):
    pi, kappa, omega, tau, edge_rates = unpack(X)
    print('pi:', pi)
    print('kappa:', kappa)
    print('omega:', omega)
    print('tau:', tau)
    print('edge rates:')
    for rate in edge_rates:
        print(rate)
    print()


def get_codon_distn_and_root_prior(codon_residue_pairs, pi):
    codon_weights = np.zeros(61)
    for i, (codon, r) in enumerate(codon_residue_pairs):
        codon_weights[i] = np.prod([pi['ACGT'.index(x)] for x in codon])
    codon_distn = codon_weights / codon_weights.sum()
    root_prior = dict(
            states = [[i, i] for i in range(61)],
            probabilities = codon_distn.tolist())
    return codon_distn, root_prior


def objective(scene, codon_residue_pairs, X):

    # For more verbosity print the parameter values.
    print_packed(X)

    # Define the log likelihood request and copy the scene.
    log_likelihood_request = dict(property = "SNNLOGL")
    scene = copy.deepcopy(scene)

    # Decode the parameter values from their unbounded
    # and unconstrained representation.
    pi, kappa, omega, tau, edge_rates = unpack(X)

    # Compute the codon distribution and root prior.
    codon_distn, root_prior = get_codon_distn_and_root_prior(
            codon_residue_pairs, pi)

    # Define the stochastic process.
    defn = get_geneconv_process_definition(
            pi, kappa, omega, tau, codon_distn, codon_residue_pairs)

    # Set the edge rates and the process definition.
    scene['root_prior'] = root_prior
    scene['tree']['edge_rate_scaling_factors'] = edge_rates
    scene['process_definitions'] = [defn]
    j_in = dict(
            scene = scene,
            requests = [log_likelihood_request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    ll = j_out['responses'][0]

    print(j_out)
    print(ll)

    return -ll


def objective_and_gradient(scene, codon_residue_pairs, X):

    # For more verbosity print the parameter values.
    print_packed(X)

    # Initialize the delta for numerical gradient calculation.
    delta = 1e-6

    #nsites = len(scene['observed_data']['iid_observations'])

    # Define the log likelihood request and the derivatives request,
    # and copy the scene.
    log_likelihood_request = dict(property = "SNNLOGL")
    edge_gradient_request = dict(property = "SDNDERI")
    scene = copy.deepcopy(scene)

    # Decode the parameter values from their unbounded
    # and unconstrained representation.
    # Compute the codon distribution and root prior.
    # Define the stochastic process.
    pi, kappa, omega, tau, edge_rates = unpack(X)
    codon_distn, root_prior = get_codon_distn_and_root_prior(
            codon_residue_pairs, pi)
    defn = get_geneconv_process_definition(
            pi, kappa, omega, tau, codon_distn, codon_residue_pairs)

    # Set the edge rates and the process definition.
    scene['root_prior'] = root_prior
    scene['tree']['edge_rate_scaling_factors'] = edge_rates
    scene['process_definitions'] = [defn]
    j_in = dict(
            scene = scene,
            requests = [
                log_likelihood_request,
                edge_gradient_request])

    print('computing log likelihood and analytical derivatives...')
    j_out = jsonctmctree.interface.process_json_in(j_in)
    neg_log_likelihood = -j_out['responses'][0]
    edge_gradient = j_out['responses'][1]

    print(j_out)
    print(neg_log_likelihood)

    # For each non-edge-specific parameter, compute a numerical derivative.
    # This is 3 for the degrees of freedom of a nucleotide distribution,
    # plus 1 each for kappa, omega, and tau.
    k = 6
    extra_derivatives = []
    print('computing finite differences for tree-wide parameters...')
    for i in range(k):
        X2 = X.copy()
        X2[i] += delta
        pi, kappa, omega, tau, edge_rates = unpack(X2)
        codon_distn, root_prior = get_codon_distn_and_root_prior(
                codon_residue_pairs, pi)
        defn = get_geneconv_process_definition(
                pi, kappa, omega, tau, codon_distn, codon_residue_pairs)
        scene['root_prior'] = root_prior
        scene['tree']['edge_rate_scaling_factors'] = edge_rates
        scene['process_definitions'] = [defn]
        j_in = dict(
                scene = scene,
                requests = [log_likelihood_request])
        j_out = jsonctmctree.interface.process_json_in(j_in)
        neg_ll = -j_out['responses'][0]
        deriv = (neg_ll - neg_log_likelihood) / delta
        extra_derivatives.append(deriv)
        print(neg_ll, deriv)
    print()

    neg_ll_gradient = np.array(extra_derivatives + [-x for x in edge_gradient])
    return neg_log_likelihood, neg_ll_gradient


def initialization_a():
    # This is for a questionable log likelihood evaluation.
    # The log likelihood should be about -1513.003.
    # This has been independently checked a few ways,
    # including with codeml (this is possible because tau=0).

    # Hard-coded ACGT nucleotide mutational distribution.
    pi =  np.array([
        0.32427103989856332,
        0.18666711777554265,
        0.20116040714181568,
        0.28790143518407829])

    # Other hard-coded parameter values.
    kappa = 5.8695382027250913
    omega = 0.087135949678171815
    tau = 0

    # Hard-code the paralogs.
    suffix_length = 7
    paralog_to_index = {
            'YML026C' : 0,
            'YDR450W' : 1}

    # Define the filenames.
    fasta_filename = 'YML026C_YDR450W_input.fasta'
    newick_filename = 'collapsed.tree.newick'

    return (
            pi,
            kappa,
            omega,
            tau,
            suffix_length,
            paralog_to_index,
            fasta_filename,
            newick_filename,
            )


def initialization_b():
    # This is for a questionable maximum likelihood estimation.

    # Initialize some parameter values.
    pi = np.ones(4) / 4
    kappa = 2.0
    omega = 0.2
    tau = 1.0

    # Hard-code the paralogs.
    suffix_length = 7
    paralog_to_index = {
            'YLR284C' : 0,
            'YOR180C' : 1}

    # Define the filenames.
    fasta_filename = 'YLR284C_YOR180C_input.fasta'
    newick_filename = 'collapsed.tree.newick'

    return (
            pi,
            kappa,
            omega,
            tau,
            suffix_length,
            paralog_to_index,
            fasta_filename,
            newick_filename,
            )


def main():

    print('initializing...')

    # Initialize some values for one of the analyses.
    (
            pi,
            kappa,
            omega,
            tau,
            suffix_length,
            paralog_to_index,
            fasta_filename,
            newick_filename,
            ) = initialization_b()
    use_empirical_pi = True
    use_uninformative_edge_rates = True
    #use_empirical_pi = False
    #use_uninformative_edge_rates = False


    print('reading the tree...')

    # Read the tree, including leaf names and branch lengths.
    # Replace the branch lengths with some uninformative guesses.
    with open(newick_filename) as fin:
        lines = fin.readlines()
    edges, edge_rates, name_to_node = read_newick(StringIO(lines[-1]))
    edge_count = len(edges)
    node_count = edge_count + 1
    if use_uninformative_edge_rates:
        edge_rates = [0.1] * edge_count
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0] * edge_count)

    print('reading the genetic code...')

    # Define the genetic code.
    codon_residue_pairs = []
    for line in _code.splitlines():
        line = line.strip()
        if line:
            row = line.upper().split()
            idx_string, residue, codon = row
            if residue != 'STOP':
                codon_residue_pairs.append((codon, residue))
    nstates = 61
    assert_equal(len(codon_residue_pairs), nstates)
    codon_to_state = {c : i for i, (c, r) in enumerate(codon_residue_pairs)}

    print('reading the fasta file...')

    # Read the fasta file.
    # At the same time, compute the empirical nucleotide distribution
    # without regard to the underlying tree.
    acgt_counts = np.zeros(4)
    observable_nodes = []
    sequences = []
    variables = []
    with open(fasta_filename) as fin:
        lines = [line.strip() for line in fin]
        lines = [line for line in lines if line]
    for name_line, sequence_line in grouper(2, lines):
        for c in sequence_line.upper():
            acgt_counts['ACGT'.index(c)] += 1
        assert_(name_line.startswith('>'))
        suffix = name_line[-suffix_length:]
        name = name_line[1:-suffix_length]
        paralog_idx = paralog_to_index[suffix]
        sequence = []
        for triple in grouper(3, sequence_line):
            codon = ''.join(triple)
            state = codon_to_state[codon]
            sequence.append(state)
        variables.append(paralog_idx)
        observable_nodes.append(name_to_node[name])
        sequences.append(sequence)

    if use_empirical_pi:
        print('computing the empirical nucleotide distribution...')

        # Define the empirical nucleotide distribution.
        pi = acgt_counts / acgt_counts.sum()
        print('empirical nucleotide distribution:', pi)

    print('defining the distribution over codons...')

    # Define the distribution over codons.
    codon_weights = np.zeros(nstates)
    for i, (codon, r) in enumerate(codon_residue_pairs):
        codon_weights[i] = np.prod([pi['ACGT'.index(x)] for x in codon])
    codon_distribution = codon_weights / codon_weights.sum()
    root_prior = dict(
            states = [[i, i] for i in range(nstates)],
            probabilities = codon_distribution.tolist())

    print('defining the codon gene conversion process...')

    # Define the process.
    process_definition = get_geneconv_process_definition(
            pi, kappa, omega, tau, codon_distribution, codon_residue_pairs)

    print('defining the observed data...')

    # Define the observed data.
    columns = zip(*sequences)
    observed_data = dict(
            nodes = observable_nodes,
            variables = variables,
            iid_observations = [list(column) for column in columns])

    print('assembling the scene...')
    
    # Assemble the scene.
    scene = dict(
            node_count = node_count,
            process_count = 1,
            state_space_shape = [nstates, nstates],
            tree = tree,
            root_prior = root_prior,
            process_definitions = [process_definition],
            observed_data = observed_data)

    print('computing the log likelihood...')

    # Ask for the log likelihood, summed over sites.
    request = dict(property = 'SNNLOGL')
    j_in = dict(
            scene = scene,
            requests = [request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    print(j_out)

    # TODO ... do a few EM-like iterations
    # to improve the estimates of the edge rate scaling factors.
    # These iterations can be quite generic.
    """
    0) As a function of only the rate matrix,
       compute the exit rate from each state.
    1) For each edge, compute the posterior expected number of transitions,
       of all types.
    2) For each edge, compute the posterior expected exit rate.
    """

    print('computing the maximum likelihood estimates...')

    # Compute the maximum likelihood estimates.
    jac = True
    f = partial(objective_and_gradient, scene, codon_residue_pairs)
    x0 = pack(pi, kappa, omega, tau, edge_rates)
    print('beginning the search...')
    result = scipy.optimize.minimize(f, x0, jac=jac, method='L-BFGS-B')

    print(result)
    xopt = result.x

    print_packed(xopt)


main()
