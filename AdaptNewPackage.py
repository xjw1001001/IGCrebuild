"""
This is a lazy code trying to adapt Alex Grinffing's new Likelihood & derivatives
calculation package without changing almost anything from his code:
https://github.com/argriffing/ctmcaas/blob/master/adv-log-likelihoods/mg-geneconv-mle.py
"""
import json
import sys
import os
import numpy as np
import scipy.optimize
from operator import mul
import networkx as nx
import pylab
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg
import argparse
import functools
from itertools import izip_longest
import cProfile

from Bio import Phylo
from Bio import SeqIO
from cStringIO import StringIO
from Bio.Phylo import PhyloXML

pwd = os.getcwd()
sys.path.insert(1, pwd + '/AlexPackage/ctmcaas/adv-log-likelihoods')
sys.path.insert(1, pwd + '/AlexPackage/jsonctmctree/jsonctmctree')
sys.path.insert(1, pwd + '/AlexPackage/npctmctree/npctmctree')

import mg94geneconv
import mle_geneconv_common

from mg94geneconv import (
    MG94_GENECONV_Abstract,
    MG94_GENECONV_Concrete)

from mg_geneconv_common import (
        ad_hoc_fasta_reader,
        get_tree_info_with_outgroup)

# Official Python itertools recipe.
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def nts_to_codons(sequence):
    return [''.join(nts) for nts in grouper(sequence, 3)]

def get_tree_info( newicktree, seqfile, paralog, SpecAfterDupli_node = 'N1'):
    tree = Phylo.read(newicktree,"newick")
    #set node number for nonterminal nodes and specify root node
    numNode=0
    for clade in tree.get_nonterminals():
        clade.name = 'N'+str(numNode)
        numNode = numNode+1
    # Promote the basic tree to PhyloXML
    tree_phy = tree.as_phyloxml(rooted = 'True')
    
    # Make a lookup table for sequences
    fastaseqs = SeqIO.parse(open(seqfile,"rU"),'fasta')
    lookup = dict((rec.id,str(rec.seq)) for rec in fastaseqs)

    for clade in tree_phy.get_terminals():
        key = clade.name
        paralog1 = paralog[0]
        paralog2 = paralog[1]
        clade.sequences.append(lookup[key+paralog1])
        if lookup.has_key(key+paralog2):
            clade.sequences.append(lookup[key+paralog2])
        else:
            print("Clade",key,"doesn't have "+ paralog2 +" gene")

    tree_nx = Phylo.to_networkx(tree_phy)
    triples = ((u.name, v.name, d['weight']) for (u,v,d) in tree_nx.edges(data=True))
    
    T = nx.DiGraph()
    for va, vb, blen in triples:
        edge = (va, vb)
        T.add_edge(*edge)    

    
    leaves = set(v for v, degree in T.degree().items() if degree == 1)
    Outgroup = list(set(leaves).difference(nx.descendants(T, SpecAfterDupli_node)))
    print 'Outgroups are :', Outgroup

    return T, tree_phy.root.name, leaves, Outgroup

def get_data(fasta_file, nnsites = None, suffix_len = 7):
    with open(fasta_file) as fin:
        name_seq_pairs = ad_hoc_fasta_reader(fin)
        
    # Convert from nucleotide sequences to codon sequences.
    name_seq_pairs = [
            (name, nts_to_codons(seq)) for name, seq in name_seq_pairs]


    # Throttle the number of sites if requested.
    if nnsites is None:
        nsites = len(name_seq_pairs[0][1])
    else:
        nsites = nnsites
        name_seq_pairs = [(n, s[:nsites]) for n, s in name_seq_pairs]
    print('number of sites to be analyzed:', nsites)

    name_to_seq = dict(name_seq_pairs)

    observable_names = [name for name, seq in name_seq_pairs]
    observable_suffixes = [name[-suffix_len:] for name in observable_names]
    suffix_to_axis = {n:i for (i, n) in enumerate(list(set(observable_suffixes))) }
    observable_axes = [suffix_to_axis[s] for s in observable_suffixes]

    # Define the map from codon to observation index.
    codon_to_state = dict()
    for i, (codon, aa) in enumerate(mg94geneconv._gen_codon_aa_pairs()):
        codon_to_state[codon.upper()] = i

    iid_observations = []
    for site in range(nsites):
        observations = []
        for name in observable_names:
            observation = codon_to_state[name_to_seq[name][site]]
            observations.append(observation)
        iid_observations.append(observations)

    return iid_observations, observable_suffixes,observable_axes, observable_names

def main(args):
    iid_observations, paralog,observable_axes, observable_names = get_data(args.fasta, args.nsites, suffix_len)
    T, root, leaves, Outgroup = get_tree_info( args.newicktree, args.fasta, paralog, SpecAfterDupli_node = 'N1' )


    # Do something about the node representations.
    # They need to be numeric.
    # That is easy enough, I think they can be arbitrarily numbered.
    names = list(T)
    name_to_node = dict((n, i) for i, n in enumerate(names))
    edges = list(T.edges())
    edge_to_eidx = dict((e, i) for i, e in enumerate(edges))

    observable_nodes = [name_to_node[n[:-suffix_len]] for n in observable_names]

    tree_row = [name_to_node[na] for na, nb in edges]
    tree_col = [name_to_node[nb] for na, nb in edges]
    tree_process = [0 if e == ('N0', Outgroup[0]) else 1 for e in edges]

    
    # define the process associated with an initial guess
    M = MG94_GENECONV_Abstract()
    guess = M.instantiate()
    guess.set_kappa(2.0)
    guess.set_omega(0.5)
    guess.set_tau(2.0)
    guess.set_nt_probs([0.30, 0.25, 0.20, 0.25])

    args.clock = True
    if args.clock:
        edge_rates = [0.5] * (len(edges) / 2 + 1)
    else:
        # define the initial edge rate guesses
        edge_rates = [0.1] * len(edges)

    # create the initial unconstrained guess vector
    x_process = guess.get_x()
    x_rates = np.log(np.array(edge_rates))
    x = np.concatenate((x_process, x_rates))
    # define the source of the log likelihood evaluation
    fn = mle_geneconv_common.eval_ll_v3module

    # define the function to minimize
    if args.clock:
        f = functools.partial(
            #mle_geneconv_common.objective,
            Clock_wrap,
            M,
            fn,
            tree_row, tree_col, tree_process,
            observable_nodes, observable_axes, iid_observations,
            edges)

        # Bound for process rates
        bnds = [(None, None)] * 4        
        bnds.extend([(None, None)] * 3)
        # Bound for Lr rates
        bnds.extend([(None, None)])
        bnds.extend([(None, 0)] * (len(edge_rates) - 1))

        result = scipy.optimize.minimize(
            f, x, jac=True, method='L-BFGS-B', bounds = bnds)  #, callback=cb_clock
    else:
        f = functools.partial(
            #mle_geneconv_common.objective,
            mle_geneconv_common.objective_and_gradient,
            M,
            fn,
            tree_row, tree_col, tree_process,
            observable_nodes, observable_axes, iid_observations,
            edges)

    # do the search
        result = scipy.optimize.minimize(
            f, x, jac=True, method='L-BFGS-B', callback=cb)
    
    print(result)
    SaveResult(args.output, result, args.clock)

def SaveResult(output, result, clock):
    with open(output + '_Clock_' + clock + '.txt', 'w+') as f:
        f.write('Clock is ' + clock + '\n')
        for i in result:
            f.write(i + ':  ' + str(result[i]) + '\n')
        
def drawtree(T):
    nx.draw_networkx(T)
    pylab.show()

# print intermediate results in the callback
def cb(x):
    k = len(edges)
    x_process, x_edge = x[:-k], x[-k:]
    edge_rates = np.exp(x_edge)
    m = M.instantiate(x_process)
    print('kappa:', m.kappa)
    print('omega:', m.omega)
    print('tau:', m.tau)
    print('nt_probs:', m.nt_probs)
    print('edge rates:')
    for edge, rate in zip(edges, edge_rates):
        print(edge, ':', rate)
    print()
    
def cb_clock(x):
    l = len(edges) / 2 + 1
    x_process, x_edge = x[:-l], x[-l:]
    edge_rates = np.exp(x_edge)
    m = M.instantiate(x_process)
    print('kappa:', m.kappa)
    print('omega:', m.omega)
    print('tau:', m.tau)
    print('nt_probs:', m.nt_probs)
    print('Lr rates:', edge_rates)
    print()

def Clock_wrap(
        abstract_model,
        fn,
        tree_row, tree_col, tree_process,
        observable_nodes, observable_axes, iid_observations,
        edges,
        x_clock): # x_clock = [x_process, log_Lr], Lr is and array of L, r0, r1, ...

    """
    This is trying to wrap up clock constraint of branch length using Alex package
    """

    nEdge = len(edges)  # number of edges
    l = nEdge / 2 + 1   # number of leaves
    k = l - 1  # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

    x_process, Lr = x_clock[:-l], np.exp(x_clock[-l:])
    x_edge_clock = []

    for edge in edges:
        if edge[0] == 'N0':  # here I abondoned root node
            if str.isdigit(edge[1][1:]):  # (N0, N1) branch
                x_edge_clock.append( Lr[0] * Lr[1] * (1 - Lr[2]) )  # L * r0 * (1 - r1) 
            else:
                x_edge_clock.append( Lr[0] * (2 - Lr[1]) )  # L * (2 - r0)

        else:
            tmp_k = int(edge[0][1:])
            if str.isdigit( edge[1][1:] ): # ( N_temp_k, N_temp_k+1 ) branch
                x_edge_clock.append( Lr[0] * reduce( mul, Lr[1 : (tmp_k + 1)], 1)  * (1 - Lr[tmp_k + 2]) )
            else:  # ( N_temp_k, leaf ) branch
                x_edge_clock.append( Lr[0] * reduce( mul, Lr[1 : (tmp_k + 1) ], 1) )


    x_rates_clock = np.log(np.array(x_edge_clock))
    x = np.concatenate([x_process, x_rates_clock])

    print 'Lr rates: ', Lr
    print 'Process rates: ', np.exp(x_process)

    f, g = mle_geneconv_common.objective_and_gradient(
        abstract_model,
        fn,
        tree_row, tree_col, tree_process,
        observable_nodes, observable_axes, iid_observations,
        edges,
        x)

    # Now need to calculate the derivates
    other_derives, edge_derives = g[:-nEdge], g[-nEdge:]
    edge_to_derives = {edges[i] : edge_derives[i] for i in range(len(edges))}

    leaf_branch = [edge for edge in edges if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
    out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
    internal_branch = [x for x in edges if not x in leaf_branch]
    assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

    leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

    Lr_derives = []  # used to store derivatives for the clock parameters L, r0, r1, ...
    Lr_derives.append(sum(edge_derives))  # dLL/dL = sum(all derives)
    Lr_derives.append(edge_to_derives[out_group_branch] * 2 / (Lr[1] - 2)
                      + sum(edge_derives))

    for i in range(2, len(Lr)):
        Lr_derives.append( edge_to_derives[('N' + str(i - 2), 'N' + str(i - 1))] / (Lr[i] - 1)
                           + sum([edge_to_derives[internal_branch[j]] for j in range(i - 1, len(internal_branch))])  # only sum over nodes decendent from node i-1
                           + sum([edge_to_derives[leaf_branch[j]] for j in range(i - 1, len(leaf_branch))]))  # only sum over nodes decendent from node i-1

    #TODO: Need to change the two sums if using general tree

    g_clock = np.concatenate( (np.array(other_derives), np.array(Lr_derives)))

    return f, g_clock

##    f = mle_geneconv_common.objective(
##        abstract_model,
##        fn,
##        tree_row, tree_col, tree_process,
##        observable_nodes, observable_axes, iid_observations,
##        edges,
##        x)
##    
##    return f


if __name__ == '__main__':
    fasta_file = './PairsAlignemt/YDR502C_YLR180W/YDR502C_YLR180W_input.fasta'
    newicktree = './PairsAlignemt/YeastTree.newick'
    Output_file = './YDR502C_YLR180W/YDR502C_YLR180W_Codon'
    suffix_len = 7
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--ll_url')
    parser.add_argument('--newicktree', default = newicktree)
    parser.add_argument('--nsites', type=int,
            help='upper limit on the number of sites to be used')
    parser.add_argument('--fasta', default = fasta_file,
            help='fasta file with paralog alignment of selected paralog pairs')
    parser.add_argument('--clock', default = False,
            help='control of molecular clock constraint')
    parser.add_argument('--output', default = Output_file,
            help='Output file')
    #cProfile.run('main(parser.parse_args())')
    main(parser.parse_args())

    


    

    
    
    
