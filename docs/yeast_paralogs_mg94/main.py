from __future__ import print_function, division

import itertools
from StringIO import StringIO

import numpy as np
from numpy.testing import assert_equal, assert_

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

def read_newick(fin):
    # use dendropy to read this newick file
    t = dendropy.Tree(stream=fin, schema='newick')
    nodes = list(t.preorder_node_iter())

    # node index lookup
    id_to_idx = {id(n) : i for i, n in enumerate(nodes)}

    # build the networkx tree
    edges = []
    edge_rates = []
    for dendro_edge in t.preorder_edge_iter():
        if dendro_edge.tail_node and dendro_edge.head_node:
            na = id_to_idx[id(dendro_edge.tail_node)]
            nb = id_to_idx[id(dendro_edge.head_node)]
            edges.append((na, nb))
            edge_rates.append(dendro_edge.length)

    # Map leaf names to nodes.
    name_to_node = {str(n.taxon) : id_to_idx[id(n)] for n in t.leaf_nodes()}

    return edges, edge_rates, name_to_node

def codon_to_triple(codon):
    return ['ACGT'.index(c) for c in codon]

def gen_mg94_structure(codon_residue_pairs):
    # Yield (i, j, ts, tv, non, syn, nt).
    transitions = (('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C'))
    for i, (ca, ra) in enumerate(codon_residue_pairs):
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
                    yield i, j, ts, tv, non, syn, nt


def main():

    # Hard-coded ACGT nucleotide mutational distribution.
    pi =  np.array([
        0.32427103989856332,
        0.18666711777554265,
        0.20116040714181568,
        0.28790143518407829])

    # Other hard-coded parameter values.
    kappa = 5.8695382027250913
    omega = 0.087135949678171815

    # Read the tree, including leaf names and branch lengths.
    with open('yeast.tree.newick') as fin:
        lines = fin.readlines()
    edges, edge_rates, name_to_node = read_newick(StringIO(lines[-1]))
    edge_count = len(edges)
    node_count = edge_count + 1
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0] * edge_count)

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

    # Define the distribution over codons.
    codon_weights = np.zeros(nstates)
    for i, (codon, r) in enumerate(codon_residue_pairs):
        codon_weights[i] = np.prod([pi['ACGT'.index(x)] for x in codon])
    codon_distribution = codon_weights / codon_weights.sum()
    root_prior = dict(
            states = [[i] for i in range(nstates)],
            probabilities = codon_distribution.tolist())

    # Define the MG94 process.
    row_states = []
    column_states = []
    exit_rates = np.zeros(nstates)
    rates = []
    for i, j, ts, tv, non, syn, nt in gen_mg94_structure(codon_residue_pairs):
        row_states.append([i])
        column_states.append([j])
        rate = (kappa * ts + tv) * (omega * non + syn) * pi[nt]
        exit_rates[i] += rate
        rates.append(rate)
    expected_rate = np.dot(codon_distribution, exit_rates)
    transition_rates = [r / expected_rate for r in rates]
    process_definition = dict(
            row_states = row_states,
            column_states = column_states,
            transition_rates = transition_rates)

    # Read the fasta file.
    observable_nodes = []
    sequences = []
    with open('YML026C_YDR450W_input.fasta') as fin:
        lines = [line.strip() for line in fin]
        lines = [line for line in lines if line]
    for name_line, sequence_line in grouper(2, lines):
        assert_(name_line.startswith('>'))
        name = name_line[1:]
        sequence = []
        for triple in grouper(3, sequence_line):
            codon = ''.join(triple)
            state = codon_to_state[codon]
            sequence.append(state)
        observable_nodes.append(name_to_node[name])
        sequences.append(sequence)

    # Define the observed data.
    columns = zip(*sequences)
    observed_data = dict(
            nodes = observable_nodes,
            variables = [0] * len(observable_nodes),
            iid_observations = [list(column) for column in columns])
    
    # Assemble the scene.
    scene = dict(
            node_count = node_count,
            process_count = 1,
            state_space_shape = [nstates],
            tree = tree,
            root_prior = root_prior,
            process_definitions = [process_definition],
            observed_data = observed_data)

    # Ask for the log likelihood, summed over sites.
    request = dict(property = 'SNNLOGL')

    j_in = dict(
            scene = scene,
            requests = [request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    print(j_out)


main()
