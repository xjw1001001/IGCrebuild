# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonBased2Repeats class
from operator import mul
from itertools import product
import os, sys

import numpy as np
import networkx as nx
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg

from Bio import Phylo
from Bio import SeqIO


pwd = os.getcwd()
sys.path.insert(1, pwd + '/AlexPackage/jsonctmctree/jsonctmctree')
import jsonctmctree.ll

class CodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, nnsites = None):
        self.newicktree = tree_newick  # newick tree file loc
        self.seqloc = alignment # multiple sequence alignment, now need to remove gap before-hand
        self.paralog = paralog  # paralog list
        self.nsites = nnsites  # number of sites in the alignment used for calculation


        # Tree topology related variable
        self.SpecAfterDupli_node = 'N1'
        self.get_tree()        
         # self.tree_phy
         # self.edge_to_blen
         # self.node_to_num
         # self.tree : {tree_row:[], tree_col:[], tree_process:[]}
        
        bases = 'tcag'.upper()
        self.nt_to_state = {a:i for (i, a) in enumerate('ACGT')}
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

        self.codon_table = dict(zip(codons, amino_acids))
        self.codon_nonstop = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.pair_to_state = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}


        # Tip data related variable
        self.get_data()
         # self.name_to_seq (actrually name_to_codonlist)
         # self.observable_names
         # self.observable_nodes
         # self.observable_axes
         # self.iid_observations

        # Rate matrix related variable
        self.get_initial_x_process()
        self.unpack_x_process()
         # self.x_process
         # self.x_rates
         # self.x (x_process + x_rates)
         # self.pi
         # self.kappa
         # self.omega
         # self.tau
        self.processes = [self.get_MG94Sparse(), self.get_MG94Geneconv()]
        self.get_prior() 

        # guess related variable

    def get_initial_x_process(self, kappa = 1.2, omega = 0.9, tau = 1.4):
        count = np.array([0, 0, 0, 0], dtype = float) # count for A, C, G, T in all seq
        for name in self.name_to_seq:
            for i in range(4):
                count[i] += ''.join(self.name_to_seq[name]).count('ACGT'[i])
        count = count / count.sum()
        
        # x_process[] = %AG, %A, %C, kappa, omega, tau
        x_process = np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                              kappa, omega, tau])
        
        setattr(CodonGeneconv, 'x_process', x_process)
        x_rates = np.array([ self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()])
        setattr(CodonGeneconv, 'x_rates', x_rates)
        x = np.concatenate((x_process, x_rates))
        setattr(CodonGeneconv, 'x', x)
        

    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        setattr(CodonGeneconv, 'tree_phy', tree_phy)
        # I decided not to append data onto topology this time
        # Because they were stored in self.name_to_seq
        tree_nx = Phylo.to_networkx(self.tree_phy)

        triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        setattr(CodonGeneconv, 'edge_to_blen', edge_to_blen)
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        Outgroup = list(set(leaves).difference(nx.descendants(T, self.SpecAfterDupli_node)))
        node_names = list(T)
        node_to_number = {n: i for i, n in enumerate(node_names)}
        setattr(CodonGeneconv, 'node_to_num', node_to_number)
        tree_row = [self.node_to_num[na] for na, nb in edge_to_blen.keys()]
        tree_col = [self.node_to_num[nb] for na, nb in edge_to_blen.keys()]
        tree_process = [0 if e == ('N0', Outgroup[0]) else 1 for e in edge_to_blen]
        tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )
        setattr(CodonGeneconv, 'tree', tree)

        
    def unpack_x_process(self):
        # x_process[] = %AG, %A, %C, kappa, omega, tau
        assert(len(self.x_process) == 6)
        x_process = self.x_process
        pi_a = x_process[0] * x_process[1]
        pi_c = (1 - x_process[0]) * x_process[2]
        pi_g = x_process[0] * (1 - x_process[1])
        pi_t = (1 - x_process[0]) * (1 - x_process[2])
        self.pi = [pi_a, pi_c, pi_g, pi_t]
        self.kappa = x_process[3]
        self.omega = x_process[4]
        self.tau = x_process[5]

    def unpack_x_rates(self):
        x_rates = self.x_rates
        assert(len(x_rates) == len(x_rates))
        for it in len(self.edge_to_blen):
            edge = self.edge_to_blen.keys()[it]
            self.edge_to_blen[edge] = x_rates[it]
        

    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(self.name_to_seq[name]) / 3 )]
            self.name_to_seq[name] = tmp_seq
            
        
    def get_data(self):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        setattr(CodonGeneconv, 'name_to_seq', {name:str(seq_dict[name].seq) for name in seq_dict.keys()})

        # Convert from nucleotide sequences to codon sequences.
        self.nts_to_codons()
        
        # Throttle the number of sites if requested.
        if self.nsites is None:
            self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])
        else:
            for name in self.name_to_seq:
                self.name_to_seq[name] = self.name_to_seq[name][: self.nsites]

        print('number of sites to be analyzed:', self.nsites)

        observable_names = self.name_to_seq.keys()
        setattr(CodonGeneconv, 'observable_names', observable_names)
        paralog_len = [len(a) for a in self.paralog]
        assert(paralog_len[1:] == paralog_len[:-1])  # check if all paralog names have same length
        suffix_len = len(self.paralog[0])
        #observable_suffixes = [name[-suffix_len:] for name in observable_names]
        observable_suffixes = self.paralog
        suffix_to_axis = {n:i for (i, n) in enumerate(list(set(observable_suffixes))) }
        observable_nodes = [self.node_to_num[n[:-suffix_len]] for n in observable_names]
        observable_axes = [suffix_to_axis[s[-suffix_len:]] for s in observable_names]
        setattr(CodonGeneconv, 'observable_axes', observable_axes)
        setattr(CodonGeneconv, 'observable_nodes', observable_nodes)


        iid_observations = []
        for site in range(self.nsites):
            observations = []
            for name in self.observable_names:
                observation = self.codon_to_state[self.name_to_seq[name][site]]
                observations.append(observation)
            iid_observations.append(observations)
        setattr(CodonGeneconv, 'iid_observations', iid_observations)

    def get_MG94Basic(self):
        # May not need to construct this matrix as well
        # unpack Pi parameters
        # Pi = %A, %C, %G, %T
        Qbasic = np.zeros((61,61),dtype = float)

        for i, ca in enumerate(self.codon_nonstop):
            for j, cb in enumerate(self.codon_nonstop):
                if i == j:
                    continue
                Qbasic[i, j] = get_MG94BasicRate(ca, cb)

        distn = [ reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
        distn = np.array(distn) / sum(distn)
        return Qbasic, distn

    def get_prior(self):
        prior_feasible_states = [(self.codon_to_state[codon], self.codon_to_state[codon]) for codon in self.codon_nonstop]
        distn = [ reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
        distn = np.array(distn) / sum(distn)
        setattr(CodonGeneconv, 'prior_feasible_states', prior_feasible_states)
        setattr(CodonGeneconv, 'prior_distribution', distn)

    def get_MG94BasicRate(self, ca, cb):
        dif = [ii for ii in range(3) if ca[ii] != cb[ii]]
        ndiff = len(dif)
        if ndiff > 1:
            return 0
        elif ndiff == 0:
            print 'Please check your codon tables and make sure no redundancy'
            print ca, cb
            return 0
        else:
            na = ca[dif[0]]
            nb = cb[dif[0]]
            QbasicRate = self.pi['ACGT'.index(nb)]

            if self.isTransition(na, nb):
                QbasicRate *= self.kappa

            if self.isNonsynonymous(ca, cb):
                QbasicRate *= self.omega

            return QbasicRate

    def isTransition(self, na, nb):
        return (set([na, nb]) == set(['A', 'G']) or set([na, nb]) == set(['C', 'T']))

    def isNonsynonymous(self, ca, cb):
        return (self.codon_table[ca] != self.codon_table[cb])

    def get_MG94Geneconv_and_MG94(self):
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            rate_sum_geneconv = 0
            rate_sum_basic = 0
            # use ca, cb, cc to denote codon_a, codon_b, codon_c, where cc != ca, cc != cb
            ca, cb = pair
            sa = self.codon_to_state[ca]
            sb = self.codon_to_state[cb]
            if ca != cb:
                for cc in self.codon_nonstop:
                    if cc == ca or cc == cb:
                        continue
                    sc = self.codon_to_state[cc]
                    # (ca, cb) to (ca, cc)
                    Qbasic = self.get_MG94BasicRate(cb, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qbasic)
                        rate_sum_geneconv += Qbasic
                        rate_basic.append(Qbasic)
                        rate_sum_basic += Qbasic

                    # (ca, cb) to (cc, cb)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate_geneconv.append(Qbasic)
                        rate_sum_geneconv += Qbasic
                        rate_basic.append(Qbasic)
                        rate_sum_basic += Qbasic

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qbasic = self.get_MG94BasicRate(cb, ca)
                if self.isNonsynonymous(cb, ca):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                rate_geneconv.append(Qbasic + Tgeneconv)
                rate_sum_geneconv += Qbasic + Tgeneconv
                rate_basic.append(Qbasic)
                rate_sum_basic += Qbasic
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qbasic = self.get_MG94BasicRate(ca, cb)
                rate_geneconv.append(Qbasic + Tgeneconv)
                rate_sum_geneconv += Qbasic + Tgeneconv
                rate_basic.append(Qbasic)
                rate_sum_basic += Qbasic

                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate_geneconv.append(-rate_sum_geneconv)
                rate_basic.append(-rate_sum_basic)
            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qbasic)
                        rate_basic.append(Qbasic)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate_geneconv.append(Qbasic)
                        rate_basic.append(Qbasic)
                        rate_sum_geneconv += 2 * Qbasic
                        rate_sum_basic += 2 * Qbasic
                        
                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate_geneconv.append(-rate_sum_geneconv)
                rate_basic.append(-rate_sum_basic)

        
    def loglikelihood_and_gradient(self):
        '''
        Edited from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        delta = 1e-8

        k = len(self.edge_to_blen)
        self.x_process, self.x_edge = self.x[:-k], self.x[-k:]

        # prepare some extra parameters for the json interface
        requested_derivatives = list(range(k))
        site_weights = np.ones(self.nsites)

        # prepare the input for the json interface
        data = dict(
            site_weights = site_weights,
            requested_derivatives = requested_derivatives,
            node_count = len(self.edge_to_blen) + 1,
            state_space_shape = [61, 61],
            process_count = len(self.processes),
            processes = self.processes,
            tree = self.tree,
            prior_feasible_states = self.prior_feasible_states,
            prior_distribution = self.prior_distribution,
            observable_nodes = self.observable_nodes,
            observable_axes = self.observable_axes,
            iid_observations = self.iid_observations
            )
        j_ll = jsonctmctree.ll.process_json_in(data)
        return j_ll
        #ll, edge_derivs = 

        


if __name__ == '__main__':
    alignment_file = '/Users/xji3/Genconv/PairsAlignemt/YDR502C_YLR180W/YDR502C_YLR180W_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    Output_file = '/Users/xji3/Genconv/YDR502C_YLR180W/YDR502C_YLR180W_Codon'
    paralog = ['YDR502C', 'YLR180W']
    test = CodonGeneconv( newicktree, alignment_file, paralog)

    test.loglikelihood_and_gradient()
            
