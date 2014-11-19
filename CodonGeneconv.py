# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonBased2Repeats class
from operator import mul
from itertools import product

import numpy as np
import networkx as nx
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg

from Bio import Phylo
from Bio import SeqIO 

class CodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, nnsites = None):
        self.newicktree = tree_newick  # newick tree file loc
        self.seqloc = alignment # multiple sequence alignment, now need to remove gap before-hand
        self.paralog = paralog  # paralog list
        self.nsites = nnsites  # number of sites in the alignment used for calculation


        # Tree topology related variable
        self.get_tree()
         # self.tree_phy
         # self.edge_to_blen
         # self.node_to_num
        
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
         # self.observable_axes
         # self.iid_observations

        # Rate matrix related variable
         # self.pi
         # self.kappa
         # self.omega
         # self.tau
        self.pi = [0.25, 0.25, 0.25, 0.25]
        self.kappa = 1.2
        self.omega = 0.9
        self.tau = 1.4
        self.processes = [self.get_MG94Sparse()]#, self.get_MG94Geneconv()]

        
        

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

        # Now create the edge_to_blen dictionary
        

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
        observable_axes = [suffix_to_axis[s] for s in observable_suffixes]
        setattr(CodonGeneconv, 'observable_axes', observable_axes)


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

    def get_MG94Geneconv(self):
        # Change idea here
        # no need to form the huge matrix but only track the non-zero entries
        #Qbasic, distn = self.get_MG94Basic()
        row = []
        col = []
        rate = []
        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            row_sum = 0
            # use ca, cb, cc to denote codon_a, codon_b, codon_c, where cc != ca, cc != cb
            ca, cb = pair
            sa = self.codon_to_state[ca]
            sb = self.codon_to_state[cb]
            for cc in self.codon_nonstop:
                if cc == ca or cc == cb:
                    continue
                sc = self.codon_to_state[cc]
                # (ca, cb) to (ca, cc)
                Qbasic = self.get_MG94BasicRate(cb, cc)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sa, sc))
                    rate.append(Qbasic)
                    row_sum += Qbasic

                # (ca, cb) to (cc, cb)
                Qbasic = self.get_MG94BasicRate(ca, cc)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sc, sb))
                    rate.append(Qbasic)
                    row_sum += Qbasic

                    
            # (ca, cb) to (ca, ca)
            row.append((sa, sb))
            col.append((sa, sa))
            Qbasic = self.get_MG94BasicRate(cb, ca)
            if self.isNonsynonymous(cb, ca):
                Tgeneconv = self.tau * self.omega
            else:
                Tgeneconv = self.tau
            rate.append(Qbasic + Tgeneconv)
            row_sum += Qbasic + Tgeneconv
            
            # (ca, cb) to (cb, cb)
            row.append((sa, sb))
            col.append((sb, sb))
            Qbasic = self.get_MG94BasicRate(ca, cb)
            rate.append(Qbasic + Tgeneconv)
            row_sum += Qbasic + Tgeneconv

            # Finally add the diagonal
            # (ca, cb) to (ca, cb)
            row.append((sa, sb))
            col.append((sa, sb))
            rate.append(-row_sum)

        process = {'row':row,
                   'col':col,
                   'rate':rate}
        return process

    def get_MG94Sparse(self):
        # Almost the same work flow with the Geneconv one because I am lazy
        row = []
        col = []
        rate = []
        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            row_sum = 0
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
                        rate.append(Qbasic)
                        row_sum += Qbasic

                    # (ca, cb) to (cc, cb)
                    Qbasic = self.get_MG94BasicRate(ca, cc)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate.append(Qbasic)
                        row_sum += Qbasic

                
                # (ca, cb) to (ca, ca)
                Qbasic = self.get_MG94BasicRate(cb, ca)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sa, sa))
                    rate.append(Qbasic)
                    row_sum += Qbasic
                
                # (ca, cb) to (cb, cb)
                Qbasic = self.get_MG94BasicRate(ca, cb)
                if Qbasic != 0:
                    row.append((sa, sb))
                    col.append((sb, sb))
                    rate.append(Qbasic)
                    row_sum += Qbasic

                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
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
                        rate.append(Qbasic)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate.append(Qbasic)

                        row_sum += 2 * Qbasic
                row.append((sa, sb))
                col.append((sa, sb))
                rate.append(-row_sum)
           

        process = {'row':row,
                   'col':col,
                   'rate':rate}
        return process
        
        


if __name__ == '__main__':
    alignment_file = '/Users/xji3/Genconv/PairsAlignemt/YDR502C_YLR180W/YDR502C_YLR180W_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    Output_file = '/Users/xji3/Genconv/YDR502C_YLR180W/YDR502C_YLR180W_Codon'
    paralog = ['YDR502C', 'YLR180W']
    test = CodonGeneconv( newicktree, alignment_file, paralog)
            
