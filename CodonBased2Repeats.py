"""
Codon based Gene Conversion rate substitution model based on Alex Griffing's npmctree package
"""

import sys

#from geneconv import *
#from __future__ import division, print_function, absolute_import

from functools import partial
from itertools import product
import argparse

import numpy as np
import networkx as nx
import pylab
import scipy
import scipy.optimize
import scipy.sparse
import scipy.sparse.linalg
#import scipy.linalg
from numpy.testing import assert_allclose

import npmctree
from npmctree.sampling import sample_histories

from npmctree.sampling import sample_histories,sample_history
from npmctree.dynamic_lmap_lhood import get_iid_lhoods, get_lhood, get_edge_to_distn2d

from Bio import Phylo
from Bio import SeqIO
from cStringIO import StringIO
from Bio.Phylo import PhyloXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from copy import deepcopy
import subprocess
import operator
from operator import mul
import pickle


class Codon2RepeatsPhy:
    def __init__(self,numLeaf,blen,tree_newick,dataloc,paralog=['EDN','ECP'], root_branch = ('N0','Tamarin'),cdmodel = False,align=False, removegaps = False):
        self.nleaf = numLeaf
        self.nbranch = 2*numLeaf -2
        self.blen = blen
        self.paralog = paralog
##        if not self.nbranch == len(blen):
##            print 'please make sure the vector blen contains all leaves'
        self.newicktree = tree_newick
        self.seqloc = dataloc
        self.nsites = 0
        self.treetopo, self.root, self.edge_to_blen, self.tree_phy = self.get_tree_info(align, removegaps)
        self.models={'Jukes-Cantor':0,'K80':1, 'F81':2, 'HKY':3,'GTR':4,'Codon':6, 'Codon_corrected':7}
        self.duplication_node = 'N0' #This is for test version, need to change later
        self.SpecAfterDupli_node = 'N1'#'Gorilla'#'N1' #This is for test version, need to change later 
        self.data = self.get_data(cdmodel)
        self.err = 1e-10
        self.para = []
        self.modelnum = 0
        self.Tao = 0.0
        self.root_branch = root_branch
        
    def DataPre(self):
        try:
            cmdline = MuscleCommandline(input=self.seqloc, out="data_out.aln", clw=True)
            cmdline()
        except:
            cmdline = ['./muscle','-in',self.seqloc, '-out', "data_out.aln", '-clw']
            subprocess.check_output(cmdline)
            
        align = AlignIO.read("data_out.aln","clustal")
        #Now, remove gaps
        i=0
        while(i<align.get_alignment_length()):
            if not align[:,i].find('-') == -1:
                if i==0:
                    align = align[:,1:]
                elif i==(align.get_alignment_length()-1):
                    align = align[:,:-1]
                else:
                    align = align[:,:i]+align[:,(i+1):]
            else:
                i=i+1
        self.nsites = align.get_alignment_length()
        return align

    def PrintResult(self, output, r1):
        with open(output,'w+') as f: #'./RootedTest_HKY.txt'
            f.write('Model Number is '+str(self.modelnum)+'\t\n')
            f.write('Parameter estimates are ')
            for item in self.para:
                f.write(str(item)+'\t')
            f.write('\n')
            f.write('Estimated Tau is ' + str(self.Tao))
            f.write('\t\n')
            f.write('Estimated branch lengths are')
            for item in self.edge_to_blen:
                f.write(str(item)+'\t'+str(self.edge_to_blen[item])+'\n')
            for i in r1:
                f.write(i+':  ' +str(r1[i])+'\n')
        
     

    def get_tree_info(self,align = False, removegaps = False):
        tree = Phylo.read(self.newicktree,"newick")
        #set node number for nonterminal nodes and specify root node
        numNode=0
        for clade in tree.get_nonterminals():
            clade.name = 'N'+str(numNode)
            numNode = numNode+1
        # Promote the basic tree to PhyloXML
        tree_phy = tree.as_phyloxml(rooted = 'True')
        # Make a lookup table for sequences
        #fastaseqs = SeqIO.parse(open(self.seqloc,"rU"),'fasta')
        if align:
            lookup = dict((rec.id,str(rec.seq)) for rec in self.DataPre())
        elif removegaps:
            fastaseqs = SeqIO.parse(open(self.seqloc,"rU"),'fasta')
            lookup = dict((rec.id,str(rec.seq)) for rec in fastaseqs)
            gap_positions = []
            for k in lookup.keys():
                gap_positions.extend([i for i,x in enumerate(lookup[k]) if x=='-'])
            gap_positions = set(gap_positions)
            for k in lookup.keys():
                for index in sorted(gap_positions, reverse = True):
                     lookup[k] = lookup[k][:index]+lookup[k][index+1:]
            self.nsites = len(lookup[lookup.keys()[0]])
        else:
            fastaseqs = SeqIO.parse(open(self.seqloc,"rU"),'fasta')
            lookup = dict((rec.id,str(rec.seq)) for rec in fastaseqs)
            self.nsites = len(lookup[lookup.keys()[0]])
        
        for clade in tree_phy.get_terminals():
            key = clade.name
            paralog1 = self.paralog[0]
            paralog2 = self.paralog[1]
            clade.sequences.append(lookup[key+paralog1])
            if lookup.has_key(key+paralog2):
                clade.sequences.append(lookup[key+paralog2])
            else:
                print("Clade",key,"doesn't have"+ paralog2 +" gene")
            
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name,v.name,d['weight']) for (u,v,d) in tree_nx.edges(data=True))
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen
        for i, k in enumerate(edge_to_blen):
            edge_to_blen[k] = self.blen[i]

#        edge_to_blen = {(u.name,v.name):d['weight'] for (u,v,d) in tree_nx.edges(data=True)}
        return T, tree_phy.root.name, edge_to_blen, tree_phy

    def drawtree(self):
        nx.draw_networkx(self.treetopo)
        pylab.show()

    def getcasenum(self,SubModel):
        #By default, use Jukes-Cantor Model
        casenum=0
        if type(SubModel) == str:
            try:
                casenum = self.models[SubModel]
            except:
                print('Please make sure the model name is valid, please choose one from the following')
                print(self.models.keys())
        elif type(SubModel) == int:
            casenum = SubModel
        else:
            print('Submodel value is either string or int, please check')
        return casenum

    def get_data(self,codonmodel=False):
        d = {}
        name = self.paralog
        for clade in self.tree_phy.get_terminals():
            for i in range(len(clade.sequences)):
                if not codonmodel:
                    d[clade.name,name[i]] = clade.sequences[i]
                else:
                    cdlist = [clade.sequences[i][3*j:(3*j+3)].upper() for j in range(len(clade.sequences[i])/3)]
                    d[clade.name,name[i]] = cdlist
            if len(clade.sequences) == 1: 
                d[clade.name,name[1]] = d[clade.name,name[0]]
        
        return d
            
    def get_state_space(self,NumRepeats):
        nt_pairs=[]
        pair_to_state ={}
        for i, pair in enumerate(product('ACGT', repeat=NumRepeats)):
            nt_pairs.append(pair)
            pair_to_state[pair] = i
        return nt_pairs, pair_to_state

    def get_codon_state_space(self,NumRepeats):
        codon_pairs = []
        pair_to_state = {}
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, amino_acids))
        codon_nonstop = [a for a in codon_table.keys() if not codon_table[a]=='*']
        for i, pair in enumerate(product(codon_nonstop,repeat=NumRepeats)):
            codon_pairs.append(pair)
            pair_to_state[pair] = i
        return codon_pairs,pair_to_state

    def Q_norm(self,SubModel,para):
        casenum = self.getcasenum(SubModel)
        if casenum == 0:
            #By default, use Jukes-Cantor Model
            Qnorm = (np.ones((4,4),dtype=float)-np.identity(4,dtype=float)*4)*para[0]
            #para[]=k
            dist = 0.25*np.ones((1,4),dtype=float)
            
        elif casenum == 1:
            #K80 Kimura, 1980 model
            #para[]=k
            Qnorm = np.ones((4,4),dtype=float)-np.identity(4,dtype=float)
            ntlist='ACGT'
            for i,na in enumerate(ntlist):
                for j, nb in enumerate(ntlist):
                    if i==j:
                        continue
                    if set([na,nb])==set(['A','G']) or set([na,nb])==set(['C','T']):
                        Qnorm[i,j] = para[0]
            Qnorm = Qnorm - np.diag(Qnorm.sum(axis=1))
            dist = 0.25*np.ones((1,4),dtype=float)
            
        elif casenum == 2:
            #F81 model, Felsenstein 1981 model
            #para[]=%AG,%A,%C
            pi_a = para[0]*para[1]
            pi_c = (1-para[0])*para[2]
            pi_g = para[0]*(1-para[1])
            pi_t = (1-para[0])*(1-para[2])
            _pi = [pi_a,pi_c,pi_g,pi_t]
            Qnorm = np.zeros((4,4),dtype=float)
            for i in range(0,4):
                Qnorm[:,i]=_pi[i]
            for i in range(0,4):
                Qnorm[i,i]=0.0
            Qnorm = Qnorm - np.diag(Qnorm.sum(axis=1))
            dist = np.array(_pi)
        
        elif casenum == 3:
            #HKY model
            #para[]=%AG,%A,%C,k
            pi_a = para[0]*para[1]
            pi_c = (1-para[0])*para[2]
            pi_g = para[0]*(1-para[1])
            pi_t = (1-para[0])*(1-para[2])
            _pi = [pi_a,pi_c,pi_g,pi_t]
            Qnorm = np.zeros((4,4),dtype=float)
            for i in range(0,4):
                Qnorm[:,i]=_pi[i]
            ntlist='ACGT'
            for i,na in enumerate(ntlist):
                for j, nb in enumerate(ntlist):
                    if i==j:
                        Qnorm[i,j]=0.0
                    if set([na,nb])==set(['A','G']) or set([na,nb])==set(['C','T']):
                        Qnorm[i,j] = Qnorm[i,j]*para[3]
            Qnorm = Qnorm - np.diag(Qnorm.sum(axis=1))
            dist = np.array(_pi)

        elif casenum == 6 or casenum == 7:
            #Codon model
            #para[] = %AG, %A, %C, k, w
            pi_a = para[0]*para[1]
            pi_c = (1-para[0])*para[2]
            pi_g = para[0]*(1-para[1])
            pi_t = (1-para[0])*(1-para[2])
            _pi = [pi_a,pi_c,pi_g,pi_t]
            Qnorm = np.zeros((61,61),dtype = float)
            codon_pairs, pairs_to_state = self.get_codon_state_space(2)
            bases = 'tcag'.upper()
            codons = [a+b+c for a in bases for b in bases for c in bases]
            amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            codon_table = dict(zip(codons, amino_acids))
            codon_nonstop = [a for a in codon_table.keys() if not codon_table[a]=='*']

            for i,ca in enumerate(codon_nonstop):
                for j,cb in enumerate(codon_nonstop):
                    if i==j:
                        continue
                    dif = [ii for ii in range(3) if ca[ii]!=cb[ii]]
                    ndiff  = len(dif)
                    if ndiff > 1:
                        continue
                    elif ndiff == 0:
                        print 'Please check your codon tables and make sure no redundancy'
                    else:
                        na = ca[dif[0]]
                        nb = cb[dif[0]]
                        Qnorm[i,j] = _pi['ACGT'.index(nb)]
                        if set([na,nb])==set(['A','G']) or set([na,nb])==set(['C','T']):
                            Qnorm[i,j] = Qnorm[i,j]*para[3]
                        if not codon_table[ca] == codon_table[cb]:
                            Qnorm[i,j] = Qnorm[i,j]*para[4]
            Qnorm = Qnorm - np.diag(Qnorm.sum(axis=1))
            dist = np.zeros((len(codon_nonstop)))
            for i in range(len(codon_nonstop)):
                dist[i] = reduce(operator.mul, [_pi['ACGT'.index(nt)] for nt in codon_nonstop[i] ], 1)
            dist = dist/dist.sum()

            
        return Qnorm,dist
        

    def get_Q_and_distn(self,SubModel,para,Tao,NumDupli = 2):
        #By default, use Jukes-Cantor Model
        casenum = self.getcasenum(SubModel)
        Qnorm,dist = self.Q_norm(SubModel,para)
        expected_rate = np.dot(dist, -np.diag(Qnorm))
        if expected_rate == 0.0:
            print para, dist
        Qnorm = Qnorm / expected_rate
        
        
        if casenum < 6:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = self.get_codon_state_space(2)
            state_pairs = codon_pairs
            bases = 'tcag'.upper()
            codons = [a+b+c for a in bases for b in bases for c in bases]
            amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            codon_table = dict(zip(codons, amino_acids))
            codon_nonstop = [a for a in codon_table.keys() if not codon_table[a]=='*']
        nstates = len(state_pairs)
        n = len(state_pairs)
        Q_un = np.zeros((n,n),dtype=float)
        if NumDupli == 1:
            if casenum<6:
                for i, (s0a, s1a) in enumerate(state_pairs):
                    for j, (s0b, s1b) in enumerate(state_pairs):
                        if i==j:
                            continue
                        if s0a == s1a and s0b == s1b:
                            Q_un[i,j] = Qnorm['ACGT'.index(s0a),'ACGT'.index(s0b)]
                Q_un = Q_un - np.diag(Q_un.sum(axis=1))
                distn = [0.0]*n#np.zeros((1,n),dtype=float)
                for i, (a,b) in enumerate(state_pairs):
                    if a==b:
                        distn[i] = dist['ACGT'.index(a)]
                
                return Q_un,np.array(distn)
            elif casenum == 6 or casenum == 7:

                for i, (s0a, s1a) in enumerate(state_pairs):
                    for j, (s0b, s1b) in enumerate(state_pairs):
                        if i==j:
                            continue
                        if s0a == s1a and s0b == s1b:
                            Q_un[i,j] = Qnorm[codon_nonstop.index(s0a),codon_nonstop.index(s0b)]
                Q_un = Q_un - np.diag(Q_un.sum(axis=1))
                distn = [0.0]*n
                for i, (a,b) in enumerate(state_pairs):
                    if a==b:
                        distn[i] = dist[codon_nonstop.index(a)]
                return Q_un,np.array(distn)
                
        elif casenum<6:
            # Then it's nucleotide model
            for i, (s0a, s1a) in enumerate(nt_pairs):
                for j, (s0b, s1b) in enumerate(nt_pairs):
                    # Diagonal entries will be set later.
                    if i == j:
                        continue
                    # Only one change is allowed at a time.
                    if s0a != s0b and s1a != s1b:
                        continue
                    # Determine which paralog changes.
                    if s0a != s0b:
                        sa = s0a
                        sb = s0b
                        context = s1a
                    if s1a != s1b:
                        sa = s1a
                        sb = s1b
                        context = s0a
                    # Set the rate according to the kind of change.
                    if context == sb:
                        rate = Tao
                    else:
                        rate = 0.0
                    Q_un[i, j] = rate+Qnorm['ACGT'.index(sa),'ACGT'.index(sb)]
            Q_un = Q_un - np.diag(Q_un.sum(axis=1))
            distn = [0.0]*n
            for i, (a,b) in enumerate(state_pairs):
                if a==b:
                    distn[i] = dist['ACGT'.index(a)]

            
#not normalized
            return Q_un,distn

        elif casenum == 6:
            for i, (s0a, s1a) in enumerate(state_pairs):
                for j, (s0b, s1b) in enumerate(state_pairs):
                    # Diagonal entries will be set later.
                    if i == j:
                        continue
                    # Only one change is allowed at a time.
                    if s0a != s0b and s1a != s1b:
                        continue
                    # Determine which paralog changes.
                    if s0a != s0b:
                        sa = s0a
                        sb = s0b
                        context = s1a
                    if s1a != s1b:
                        sa = s1a
                        sb = s1b
                        context = s0a
                    # Set the rate according to the kind of change.
                    if context == sb:
                        rate = Tao
                    else:
                        rate = 0.0
                    Q_un[i, j] = rate+Qnorm[codon_nonstop.index(sa),codon_nonstop.index(sb)]
            Q_un = Q_un - np.diag(Q_un.sum(axis=1))
            distn = [0.0]*n
            for i, (a,b) in enumerate(state_pairs):
                if a==b:
                    distn[i] = dist[codon_nonstop.index(a)]
            return Q_un,np.array(distn)

        elif casenum == 7:
            bases = 'tcag'.upper()
            codons = [a+b+c for a in bases for b in bases for c in bases]
            amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            codon_table = dict(zip(codons, amino_acids))
            for i, (s0a, s1a) in enumerate(state_pairs):
                for j, (s0b, s1b) in enumerate(state_pairs):
                    # Diagonal entries will be set later.
                    if i == j:
                        continue
                    # Only one change is allowed at a time.
                    if s0a != s0b and s1a != s1b:
                        continue
                    # Determine which paralog changes.
                    if s0a != s0b:
                        sa = s0a
                        sb = s0b
                        context = s1a
                    if s1a != s1b:
                        sa = s1a
                        sb = s1b
                        context = s0a
                    # Set the rate according to the kind of change.
                    if context == sb:
                        if codon_table[sa]==codon_table[sb]:
                            rate = Tao*para[4]
                        else:
                            rate = Tao
                    else:
                        rate = 0.0
                    Q_un[i, j] = rate+Qnorm[codon_nonstop.index(sa),codon_nonstop.index(sb)]
            Q_un = Q_un - np.diag(Q_un.sum(axis=1))
            distn = [0.0]*n
            for i, (a,b) in enumerate(state_pairs):
                if a==b:
                    distn[i] = dist[codon_nonstop.index(a)]
            return Q_un,np.array(distn)
                    
        else:
            print('Please Check get_Q_and_distn function and make sure the case is considered')
            return Qnorm, distn
        
    def get_P(self,edge_to_blen_infer,SubModel,para,Tao):
        edge_to_P = {}
        try:
            all_dupli_nodes = nx.descendants(self.treetopo,self.SpecAfterDupli_node)
        except:
            all_dupli_nodes = set([''])
            #print 'Warning : No Duplication event specified'
        all_dupli_nodes.add(self.SpecAfterDupli_node)
        for edge in self.treetopo.edges():
            blen = edge_to_blen_infer[edge]
            
            if edge[1] in all_dupli_nodes:
                Q,dist = self.get_Q_and_distn(SubModel,para,Tao,2)
            else:
                Q,distn = self.get_Q_and_distn(SubModel,para,Tao,1)
            P = scipy.sparse.linalg.expm(float(blen) * Q)
            edge_to_P[edge] = P
        root_distn = distn
        return edge_to_P, root_distn

    def get_expected_rate(self, Q_un):
        w,v = scipy.sparse.linalg.eigs(Q_un.T, k=1, which = 'SM')
        weights = v[:,0].real
        distn = weights / weights.sum()
        expected_rate = np.dot(distn, -np.diag(Q_un))
        return expected_rate
            

    def objective(self,SubModel, T, root,nt_pairs,constraints,edge_to_blen_infer,para,Tao):
        edge_to_P, root_distn = self.get_P(edge_to_blen_infer,SubModel,para,Tao)

##        neg_ll = 0
##        for node_to_data_fvec1d in constraints:
##            lhood = get_lhood(T, edge_to_P, root, root_distn, node_to_data_fvec1d)
##            neg_ll -= np.log(lhood)

        lhoods = get_iid_lhoods(T,edge_to_P, root, root_distn, constraints)
        neg_ll = -np.log(lhoods).sum()

        return neg_ll
    def get_data_constraints(self,nodes, pair_to_state, nsites, data):
        """
        Get sequence data constraints for the likelihood calculation.

        Also return some match/mismatch counts.

        Parameters
        ----------
        nodes : set
            set of all nodes in the tree
        pair_to_state : dict
            maps nucleotide pairs to an integer state index
        nsites : integer
            number of sites in the alignment
        data : dict
            map from (node, paralog) to nucleotide sequence string

        Returns
        -------
        per_site_constraints : sequence
            for each site, a map from
        npaired_yes : integer
            number of site at which aligned paralog nucleotides match
        npaired_no : integer
            number of site at which aligned paralog nucleotides mismatch

        """
        # maybe do not hard code this...
        paralogs = self.paralog
        # assume we only have sequence data at the leaves...
        leaves = set(taxon for taxon, paralog in data)
        leaf_to_paralogs = {}
        for taxon,paralog in data:
            try:
                leaf_to_paralogs[taxon].append(data[(taxon,paralog)])
            except:
                leaf_to_paralogs[taxon]=[ data[(taxon,paralog)] ]
        nstates = len(pair_to_state)
        non_leaves = set(nodes) - set(leaf_to_paralogs)
        per_site_constraints = []
        npaired_yes = 0
        npaired_no = 0
        for site in range(nsites):
            node_to_data_fset = {}
            for node in non_leaves:
                node_to_data_fset[node] = np.ones(nstates, dtype=bool)
            for node in leaves:
                try:
                    nt0 = data[(node, paralogs[0])][site]
                except:
                    print node, data.keys()
                    print site
                    print data[(node,paralogs[0])]
                nt1 = data[node, paralogs[1]][site]
                pair = (nt0, nt1)
                if nt0 == nt1:
                    npaired_yes += 1
                else:
                    npaired_no += 1
                state = pair_to_state[pair]
                fset = np.zeros(nstates, dtype=bool)
                fset[state] = True
                node_to_data_fset[node] = fset
            per_site_constraints.append(node_to_data_fset)
        return per_site_constraints, npaired_yes, npaired_no        


    def estimate(self,args,SubModel,guess,output='./result.txt',est_blen = False, unrooted = False, force_tau = False,clock = False, print_result = False):
##        self.treetopo, self.root, self.edge_to_blen, self.tree_phy = self.get_tree_info()
##        T, root, edge_to_blen, tree_phy = test.get_tree_info()
        if clock:
            self.root='root'
            new_root_branch = ('root', self.root_branch[1])
            self.treetopo.remove_edge(*self.root_branch)
            self.treetopo.add_node('root')
            self.treetopo.add_edges_from([('root',self.root_branch[0]),('root',self.root_branch[1])])
            self.edge_to_blen.pop(self.root_branch)
            self.edge_to_blen[('root',self.root_branch[0])] = self.blen[0]
            self.edge_to_blen[('root',self.root_branch[1])] = self.blen[0]
            self.root_branch = new_root_branch
            leaves = set(v for v, degree in self.treetopo.degree().items() if degree == 1) 
        casenum = self.getcasenum(SubModel)
        if casenum < 6:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = self.get_codon_state_space(2)
            state_pairs = codon_pairs
        nodes = set(self.treetopo)
        if casenum<6:
            nsts = self.nsites
        else:
            nsts = self.nsites/3

        #print len(pair_to_state),casenum
        constraint_info = self.get_data_constraints(nodes, pair_to_state, nsts, self.data)
        constraints, npaired_yes, npaired_no = constraint_info

        if est_blen:
            one_bnd = (None,None)
            if clock:
                if casenum == 0:            
                    bnds = [one_bnd]
                    bnds.extend([(None,0.0)]*(len(leaves)-1+2))
                elif casenum == 1:
                    bnds = [one_bnd]
                    bnds.extend([(None,0.0)]*(len(leaves)-1+2))
                elif casenum == 2:
                    bnds = [one_bnd]
                    bnds.extend([(None,0.0)]*(len(leaves)-1))
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None)])
                elif casenum == 3:
                    bnds = [one_bnd]
                    bnds.extend([(None,0.0)]*(len(leaves)-1))
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None),(0, None)])
                elif casenum == 6 or casenum == 7:
                    bnds = [one_bnd]
                    bnds.extend([(None,0.0)]*(len(leaves)-1))
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None),(None, None),(0, None)])

                
                
            else:
                if casenum == 0:            
                    bnds = [one_bnd]*(len(self.edge_to_blen)+2)
                elif casenum == 1:
                    bnds = [one_bnd]*(len(self.edge_to_blen)+2)
                elif casenum == 2:
                    bnds = [one_bnd]*len(self.edge_to_blen)
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None)])
                elif casenum == 3:
                    bnds = [one_bnd]*len(self.edge_to_blen)
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None),(0, None)])
                elif casenum == 6 or casenum == 7:
                    bnds = [one_bnd]*len(self.edge_to_blen)
                    bnds.extend([(None, -0.05),(None, -0.05),(None, -0.05),(None, None),(None, None),(0, None)])

            e = deepcopy(self.edge_to_blen.keys())
            if unrooted:
                e.remove(('N0','N1'))
                bnds.pop(0)

            if force_tau:
                bnds.pop(-1)

            ff = partial(self.objective_for_blen, SubModel, state_pairs, constraints,e,unrooted,force_tau,clock, print_result)
            
                     
        else:
            f = partial(self.objective, SubModel, self.treetopo, self.root, state_pairs, constraints, self.edge_to_blen)
            if casenum == 0:
                bnds = ((0, None),(0, None))                
            elif casenum == 1:
                bnds = ((0, None),(0, None))                
            elif casenum == 2:
                bnds = ((self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, None))                
            elif casenum == 3:
                bnds = [(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, None),(0, None)]
            elif casenum == 6 or casenum == 7:
                bnds = [(None, 0.0),(None, 0.0),(None, 0.0),(None, None),(None, None),(0, None)]
            ff = lambda x:f(x[0:-1],x[-1])

        try:
            leaves = set(v for v, degree in self.treetopo.degree().items() if degree == 1)
            Outgroup = set(leaves).difference(nx.descendants(self.treetopo,self.SpecAfterDupli_node))
            print 'Outgroups are :', Outgroup
        except:
            print 'Warning : No Duplication event specified'
        print 'paralog site matches :', npaired_yes
        print 'paralog site mismatches:', npaired_no
##        nsmall = 4
##        pa = npaired_yes+0.0
##        pb = npaired_no+0.0
##        phi_hat = (pa+pb)/pb*3.0-4.0
##        print 'preliminary estimate of phi:', guess,phi_hat
        print 'negative log likelihood of preliminary estimate:', ff(guess)        
        

        result = scipy.optimize.minimize(ff, x0=guess, method='L-BFGS-B', bounds=bnds)
        self.PrintResult(output,result)
        print result
        if not clock:
            try:
                self.update_blen_phylo()
            except:
                print 'Please Check Update_blen_phylo for this case'
        else:
            print 'TODO implement update phylo tree for clock model'

        
        
        return result
    
    def objective_for_blen(self,SubModel,nt_pairs,constraints,edge_to_blen_keys,unrooted,force_tau,clock, print_result,x):
        casenum = self.getcasenum(SubModel)
        est_edge_to_blen = deepcopy(self.edge_to_blen)
        if not clock:
            for i, k in enumerate(edge_to_blen_keys):
                est_edge_to_blen[k] = np.exp(x[i])
                self.edge_to_blen[k] = np.exp(x[i])
        else:
            '''
            x[] = L, r0, r1, r2,...
            '''
            for k in edge_to_blen_keys:
                xx = np.exp(x)
                if k[0] == 'root':
                    if str.isdigit(k[1][1:]):
                        est_edge_to_blen[k] = xx[0]*(1-xx[1])
                        self.edge_to_blen[k] = xx[0]*(1-xx[1])
                    else:
                        est_edge_to_blen[k] = xx[0]
                        self.edge_to_blen[k] = xx[0]
                else:
                    tmp_k = int(k[0][1:])
                    if str.isdigit(k[1][1:]):                        
                        est_edge_to_blen[k] = xx[0]*reduce(mul, xx[1:(tmp_k+1)], 1)*(1-xx[tmp_k+2])
                        self.edge_to_blen[k] = xx[0]*reduce(mul, xx[1:(tmp_k+1)], 1)*(1-xx[tmp_k+2])
                    else:
                        est_edge_to_blen[k] = xx[0]*reduce(mul, xx[1:(tmp_k+1)], 1)
                        self.edge_to_blen[k] = xx[0]*reduce(mul, xx[1:(tmp_k+1)], 1)
                        
                

        if unrooted:
            self.edge_to_blen[('N0','N1')] = 1e-6
            est_edge_to_blen[('N0','N1')] = 1e-6
            start_of_para = len(est_edge_to_blen)-1
        else:
            if clock:
                leaves = set(v for v, degree in self.treetopo.degree().items() if degree == 1)
                start_of_para = len(leaves)
            else:
                start_of_para = len(est_edge_to_blen)

        if force_tau:
            est_Tao = 0.0
            est_para = np.exp(x[start_of_para:])
        else:
            est_para = np.exp(x[start_of_para:-1])
            est_Tao = x[-1]


        self.para = deepcopy(est_para)
        self.Tao = est_Tao
        self.modelnum = casenum

        rst = self.objective(SubModel, self.treetopo, self.root, nt_pairs, constraints, est_edge_to_blen,est_para, est_Tao)

        if print_result:
            print self.edge_to_blen
            #print 'x:', x
            print 'para :', est_para, est_Tao
            print rst

        #objective(self,SubModel, T, root,nt_pairs,constraints,edge_to_blen_infer,para,Tao)
        #print 'x:', x, rst
        

        return rst


    def simulator(self,edge_to_blen_infer,SubModel,para,Tao,nsites,codonmodel = False):
        edge_to_P, root_distn = self.get_P(edge_to_blen_infer,SubModel,para,Tao)
        node_to_data_lmap = {}
        if not codonmodel:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = self.get_codon_state_space(2)
            state_pairs = codon_pairs

        nodes = set(self.treetopo)
        leaves = set(v for v, degree in self.treetopo.degree().items() if degree == 1)
        nstates = len(state_pairs)
        for node in nodes:
            node_to_data_lmap[node] = np.ones(nstates)

        node_to_pairs = dict((node,[]) for node in nodes)
        leaf_to_seqs = dict((leaf,[]) for leaf in leaves)
        for node_to_state in sample_histories(self.treetopo, edge_to_P, self.root,
                                              root_distn, node_to_data_lmap,nsites):
            for node, state in node_to_state.items():
                state_pair = state_pairs[state]
                node_to_pairs[node].append(state_pair)
        for leaf in leaves:
            pairs = node_to_pairs[leaf]
            seqs = [''.join(x).upper() for x in zip(*pairs)]
            leaf_to_seqs[leaf].append(seqs)
        return leaf_to_seqs

    def simtofasta(self,leaf_to_seqs,out_file,paralog_name):
        with open(out_file,'w+') as f:
            for leaf in leaf_to_seqs.keys():
                f.write('>'+leaf+paralog_name[0]+'\n')
                f.write(leaf_to_seqs[leaf][0][0]+'\n')
                f.write('>'+leaf+paralog_name[1]+'\n')
                f.write(leaf_to_seqs[leaf][0][1]+'\n')
        
            
        
    def onetransition(self,blen,P,seq0, codon = False):
        if not codon:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = get_codon_state_space(2)
            state_pairs = codon_pairs
        state_num = pair_to_state[seq0]
        probabilities = P[state_num,:]
        bins = np.add.accumulate(probabilities)
        return state_pairs[np.digitize(random_sample(1), bins)]

    def dataoutpaml(self,out_loc):
        try:
            Outgroup = set(self.treetopo).difference(nx.descendants(self.treetopo,self.SpecAfterDupli_node))
        except:
            Outgroup = set([])
        paralog_name = self.paralog
        leaves = set(taxon for taxon, paralog in self.data)
        with open(out_loc+'dataoutpaml.paml','w') as f:
            f.write(str(len(leaves)+len(leaves.difference(Outgroup)))+'\t'+ str(self.nsites) + '\n')
            for spe,paralog in self.data.keys():
                if not spe in Outgroup:
                    if type(self.data[(spe,paralog)]) == str:
                        f.write(spe+paralog + '    ' + self.data[(spe,paralog)] + '\n')
                    else:
                        f.write(spe+paralog +'    ' + ''.join(self.data[(spe,paralog)]) + '\n')
                elif paralog == paralog_name[0]:
                    if type(self.data[(spe,paralog)]) == str:
                        f.write(spe+paralog + '    ' + self.data[(spe,paralog)] + '\n')
                    else:
                        f.write(spe+paralog +'    ' + ''.join(self.data[(spe,paralog)]) + '\n')

    def get_total_blen(self,outgroup_nodes = []):
        #leaves = set(v for v, degree in self.treetopo.degree().items() if degree == 1)
        try:
            nodes_involve_Outgroup = set(self.treetopo).difference(nx.descendants(self.treetopo,self.SpecAfterDupli_node))
            nodes_involve_Outgroup = nodes_involve_Outgroup.difference([self.duplication_node,self.SpecAfterDupli_node])
        except:
            nodes_involve_Outgroup = set([])
        total_blen = 0.0
        before_duplication_total_blen = 0.0
        post_duplication_total_blen = 0.0
        if nodes_involve_Outgroup:
            q,d = self.get_Q_and_distn(self.modelnum,self.para, self.Tao, 2)
            expected_rate = self.get_expected_rate(q)
            branch_involve_outgoup = [v for v in self.edge_to_blen.keys() if v[1] in nodes_involve_Outgroup]
            before_duplication_total_blen += sum(self.edge_to_blen[v] for v in branch_involve_outgoup)
            post_duplication_total_blen += expected_rate*sum(self.edge_to_blen[v] for v in set(self.edge_to_blen.keys()).difference(branch_involve_outgoup))
            total_blen = before_duplication_total_blen + post_duplication_total_blen
            ratio = post_duplication_total_blen/total_blen
            return total_blen, ratio
        else:
            branch_involve_outgoup = [v for v in self.edge_to_blen.keys() if v[1] in outgroup_nodes]
            before_duplication_total_blen += sum(self.edge_to_blen[v] for v in branch_involve_outgoup)
            total_blen += sum(self.edge_to_blen[v] for v in self.edge_to_blen.keys())
            ratio = (total_blen-before_duplication_total_blen)/total_blen
            return total_blen, ratio

    def update_blen_phylo(self,outgroup_nodes = []):
        try:
            nodes_involve_Outgroup = set(self.treetopo).difference(nx.descendants(self.treetopo,self.SpecAfterDupli_node))
            nodes_involve_Outgroup = nodes_involve_Outgroup.difference([self.duplication_node,self.SpecAfterDupli_node])
        except:
            nodes_involve_Outgroup = set([])
        edge_to_blen_normalized = deepcopy(self.edge_to_blen)
        if nodes_involve_Outgroup:
            q,d = self.get_Q_and_distn(self.modelnum,self.para, self.Tao, 2)
            expected_rate = self.get_expected_rate(q)
            branch_involve_outgoup = [v for v in self.edge_to_blen.keys() if v[1] in nodes_involve_Outgroup]
            branch_involve_duplication = set(self.edge_to_blen.keys()).difference(branch_involve_outgoup)
            for branch in branch_involve_duplication:
                edge_to_blen_normalized[branch] = edge_to_blen_normalized[branch]*expected_rate/2
        for clade in self.tree_phy.get_terminals():
            node_path = self.tree_phy.get_path(clade)
            node_path[0].branch_length = edge_to_blen_normalized[('N0',node_path[0].name)]
            for i in range(len(node_path)-1):
                node_path[i+1].branch_length = edge_to_blen_normalized[(node_path[i].name,node_path[i+1].name)]

    def get_Tau_matrix(self,Normalizing_factor = 1.0,NumRepeats = 2):
        casenum = self.getcasenum(self.modelnum)
        if casenum < 6:
            nt_pairs, pair_to_state = self.get_state_space(NumRepeats)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = self.get_codon_state_space(NumRepeats)
            state_pairs = codon_pairs
            bases = 'tcag'.upper()
            codons = [a+b+c for a in bases for b in bases for c in bases]
            amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
            codon_table = dict(zip(codons, amino_acids))
            codon_nonstop = [a for a in codon_table.keys() if not codon_table[a]=='*']

        Tau_matrix = np.zeros((len(state_pairs),len(state_pairs)))

        for i, (s0a, s1a) in enumerate(state_pairs):
            for j, (s0b, s1b) in enumerate(state_pairs):
                # Diagonal entries will be set later.
                if i == j:
                    continue
                # Only one change is allowed at a time.
                if s0a != s0b and s1a != s1b:
                    continue
                # Determine which paralog changes.
                if s0a != s0b:
                    sa = s0a
                    sb = s0b
                    context = s1a
                if s1a != s1b:
                    sa = s1a
                    sb = s1b
                    context = s0a
                # Set the rate according to the kind of change.
                if context == sb:
                    rate = self.Tao*Normalizing_factor
                else:
                    rate = 0.0
                Tau_matrix[i, j] = rate

        return Tau_matrix

    def get_edge_to_M(self,edge_to_Q,edge_to_C,edge_to_blen,NumRepeats = 2):
        # Q,d = self.get_Q_and_distn(self.modelnum, self.para, self.Tao, NumRepeats)
        # assert(Q.size == C.size)
        # Using Algorithm 3 in Tataru, Hoblth 2011 BMC Bioinformatics paper
        edge_to_P, root_distn = self.get_P(self.edge_to_blen,self.modelnum,self.para,self.Tao)
        edge_to_M= {}  # M is the uper right corner of the matrix exponential
        for edge in self.edge_to_blen:
            Q = edge_to_Q[edge]
            C = edge_to_C[edge]
            A = np.vstack( (np.hstack((Q,C)),np.hstack((np.zeros(Q.shape),Q))) )
            M = scipy.sparse.linalg.expm(float(edge_to_blen[edge]) * A)[0:Q.shape[0],Q.shape[1]:]
            P = edge_to_P[edge]
            P_zero_indices = [(np.where(P==0)[0][i],np.where(P==0)[1][i]) for i in range(len(np.where(P==0)[0]))]
            M_zero_indices = [(np.where(M==0)[0][i],np.where(M==0)[1][i]) for i in range(len(np.where(M==0)[0]))]
            for (i,j) in set(P_zero_indices).intersection(set(M_zero_indices)):
                P[i,j] += self.err
            edge_to_M[edge] = np.divide(M,P)

        return edge_to_M
    
    def get_edge_to_M_Geneconv(self, NumRepeats =2):
        Q_post_duplication,d = self.get_Q_and_distn(self.modelnum, self.para, self.Tao, NumRepeats)
        expected_rate = self.get_expected_rate(Q_post_duplication)
        normalizing_factor = NumRepeats/expected_rate
        Q_post_duplication_normalized = normalizing_factor*Q_post_duplication
        
        Q_pre_duplication,d_pre = self.get_Q_and_distn(self.modelnum, self.para, self.Tao, 1)
        Tau_matrix = self.get_Tau_matrix(normalizing_factor,NumRepeats)
        
##        Q_post_modified = deepcopy(Q_post_duplication_normalized)
##        Q_pre_modified = deepcopy(Q_pre_duplication)
##        nonzerolist = Q_post_duplication_normalized.nonzero()
##        nonzeros = [(nonzerolist[0][i],nonzerolist[1][i]) for i in range(len(nonzerolist[0]))]
##        for i in range(0,Q_post_duplication_normalized.shape[0]):
##            for j in range(0,Q_post_duplication_normalized.shape[1]):
##                if not (i,j) in nonzeros:
##                    Q_post_modified[i,j] += self.err
##                    
##        nonzerolist = Q_pre_duplication.nonzero()
##        nonzeros = [(nonzerolist[0][i],nonzerolist[1][i]) for i in range(len(nonzerolist[0]))]
##        for i in range(0,Q_pre_duplication.shape[0]):
##            for j in range(0,Q_pre_duplication.shape[1]):
##                if not (i,j) in nonzeros:
##                    Q_pre_modified[i,j] += self.err

        # Q modified is just created to avoid 0/0 case
        
        C_post = Tau_matrix  # Coeff Matrix for Geneconv events, it's been normalized
        C_post_none_conv = (Q_post_duplication_normalized - Tau_matrix)  # Coeff Matrix for none Geneconv events
        np.fill_diagonal(C_post_none_conv,0.0)
        C_pre = np.zeros(Tau_matrix.shape)
        C_pre_none_conv = deepcopy(Q_pre_duplication)
        np.fill_diagonal(C_pre_none_conv,0.0)
        

        edge_to_Q = {}
        edge_to_C = {}
        edge_to_C_none = {}
        edge_to_blen = {}
        try:
            all_dupli_nodes = nx.descendants(self.treetopo,self.SpecAfterDupli_node)
        except:
            all_dupli_nodes = set([''])
            #print 'Warning : No Duplication event specified'
        all_dupli_nodes.add(self.SpecAfterDupli_node)
        for edge in self.treetopo.edges():
            if edge[1] in all_dupli_nodes:
                edge_to_Q[edge] = Q_post_duplication_normalized
                edge_to_C[edge] = C_post
                edge_to_C_none[edge] = C_post_none_conv
                edge_to_blen[edge] = self.edge_to_blen[edge]/normalizing_factor  # ? Do I need to normalize again?
            else:
                edge_to_Q[edge] = Q_pre_duplication
                edge_to_C[edge] = C_pre
                edge_to_C_none[edge] = C_pre_none_conv
                edge_to_blen[edge] = self.edge_to_blen[edge]

    

        return self.get_edge_to_M(edge_to_Q,edge_to_C,edge_to_blen,NumRepeats),self.get_edge_to_M(edge_to_Q,edge_to_C_none,edge_to_blen,NumRepeats)

    def get_edge_to_expectednumchanges(self):
        edge_to_P, root_distn = self.get_P(self.edge_to_blen,self.modelnum,self.para,self.Tao)
        casenum = self.getcasenum(self.modelnum)
        if casenum < 6:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = self.get_codon_state_space(2)
            state_pairs = codon_pairs
        nodes = set(self.treetopo)
        if casenum<6:
            nsts = self.nsites
        else:
            nsts = self.nsites/3

        #print len(pair_to_state),casenum
        constraint_info = self.get_data_constraints(nodes, pair_to_state, nsts, self.data)
        constraints, npaired_yes, npaired_no = constraint_info

        edge_to_M,edge_to_M_none = self.get_edge_to_M_Geneconv()

        expected_geneconv = {v:0.0 for v in self.edge_to_blen.keys()}
        expected_none_geneconv = {v:0.0 for v in self.edge_to_blen.keys()}

        for node_to_data_fvec1d in constraints:
            edge_to_J = get_edge_to_distn2d(self.treetopo, edge_to_P, self.root, root_distn, node_to_data_fvec1d)
            for v in self.edge_to_blen.keys():
                expected_geneconv[v] += np.multiply(edge_to_M[v],edge_to_J[v]).sum()/nsts
                expected_none_geneconv[v] += np.multiply(edge_to_M_none[v],edge_to_J[v]).sum()/nsts

        return expected_geneconv, expected_none_geneconv

##        neg_ll = 0
##        for node_to_data_fvec1d in constraints:
##            lhood = get_lhood(T, edge_to_P, root, root_distn, node_to_data_fvec1d)
##            neg_ll -= np.log(lhood)

    def save_to_file(self, file_name = './CodonBased2RepeatsClassSave.p'):
        pickle.dump(self, open(file_name,'wb+'))  # use pickle to save the class which can be easily reconstructed by pickle.load()
                    

def cleanPAML(in_file,out_file, fasta_out = False):
    cleanedline = []
    with open(in_file, 'r') as f:
        for line in f:
            if line.strip():
                cleanedline.append(line.strip('\n'))
    with open(out_file, 'w') as f:
        f.writelines('\n'.join(cleanedline))

    if fasta_out:
        sim_paml = open(out_file,'rU')
        paralog_name = ['EDN','ECP']
        with open(out_file.replace('paml','fasta'),'w') as f:
            for record in SeqIO.parse(sim_paml,'phylip'):
                for i in range(2):
                    f.write('>'+record.id+paralog_name[i]+'\n')
                    f.write(str(rec.seq)+'\n')

  

if __name__ == '__main__':
    numLeaf = 5
    blen = np.ones([2*numLeaf-2])*2
    #blen = np.array([1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 3.0, 4.0])
    #blen = np.array([0.09, 0.1427433571102972,0.23973336228125663, 0.1304889427706489, 1e-06, 0.078655341006576895])
    tree_newick = './data/input_tree.newick'
    dataloc = './data/input_data.fasta'
    simdata = 'simdata.fasta'
    tree_newick_compare_paml = './data/input_tree_compare_paml.newick'
    dataloc_compare_paml = './data/input_data_compare_paml.fasta'
    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc,align=True)
    test3 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc,align=True)
    numLeaf = 9
    blen = np.ones([2*numLeaf-2])*2
    test2 = Codon2RepeatsPhy(numLeaf,blen,tree_newick_compare_paml,dataloc_compare_paml,align=False)
    test2.SpecAfterDupli_node = ''

    numLeaf_codon = 5
    blen_codon = np.ones([2*numLeaf_codon-2])*2
    dataloc_codon = './data/codon_alignment_nucleotide_format.fasta'
    tree_newick = './data/input_tree.newick'
    test4  = Codon2RepeatsPhy(numLeaf_codon,blen_codon,tree_newick,dataloc_codon,cdmodel= True,removegaps = True)
    
    
    sim_SubModel=3
    sim_para=[0.2,0.3,0.1,1.2]
    sim_Tao=1.0
    sim_nsites = 480
    #sim_nsites = 3000
    guess = np.log([0.7,0.5,0.4,np.e**2])
    guess = np.append(guess,sim_Tao/2.0)
    print 'Now simulate data'
    #a = test.simulator(test.edge_to_blen,sim_SubModel,sim_para,sim_Tao,sim_nsites,codonmodel=True)
    #test.simtofasta(a,simdata,['EDN','ECP'])
    #cleanPAML('./mc.paml','./simdata.paml',True)
    #test2 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,simdata, cdmodel = False)
##    test.drawtree()

    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()

##    test_main(args)

##    test.estimate(args,0)
#    test2.estimate(args,sim_SubModel,guess)

    guess_w_blen = [-1.0]*(len(test.edge_to_blen))
    guess_w_blen[0:6] = [0.0,-1.0,-2.0,-3.0,-4.0,-5.0]
    guess_w_blen.extend(guess)

    #print 'Estimate actual data'
    #test.estimate(args, sim_SubModel, guess_w_blen, True)

    print 'Now estimate simulated data'
    #r1 = test.estimate(args, sim_SubModel, guess_w_blen, est_blen = True)
    #print 'Total branch lengths and proportion of post-duplication are', test.get_total_blen()
    Phylo.write(test.tree_phy,'./HKY_Geneconv.newick','newick')
    guess_w_blen = [-2.0]*(len(test2.edge_to_blen)-1)
    guess_w_blen[0:6] = [0.0,-1.0,-2.0,-3.0,-4.0,-5.0]
    guess_w_blen.extend(guess)    
    #r2 = test2.estimate(args, sim_SubModel, guess_w_blen, est_blen = True,unrooted = True)
    #Phylo.write(test2.tree_phy,'./HKY.newick','newick')
    #print 'Total branch lengths and proportion of post-duplication are', test2.get_total_blen(['Tamarin'])
    guess_w_blen = [-2.0]*(len(test.edge_to_blen))
    guess_w_blen[0:6] = [0.0,-1.0,-2.0,-3.0,-4.0,-5.0]
    guess_w_blen.extend(guess)
    guess_w_blen.pop(-1)
    #r3 = test3.estimate(args, sim_SubModel, guess_w_blen, est_blen = True,unrooted = False,force_tau = True)
    #Phylo.write(test3.tree_phy,'./HKY_Geneconv_tau0.newick','newick')
    #print 'Total branch lengths and proportion of post-duplication are', test3.get_total_blen()
    guess = np.log([0.46789,0.27614/0.46789,0.25856/(0.25856+0.27356),2.08027,0.92449])
    guess = np.append(guess,[1.5])
    guess_w_blen = np.log([0.0195+0.0717,0.0556+0.026, 0.0136+0.0189, 0.0764+0.1361, 0.0826+0.2192, 0.0194+0.0191, 0.2188+0.1285, 0.3669])
    #These are the estimates from PAML
    guess = np.log([0.49625812,  0.58854282 , 0.48860353,  2.10255791,  1.14657581])
    guess = np.append(guess,[0.55095646406])
    guess_w_blen = np.log([0.031798041175813183, 0.034203400799371592, 0.014099750421849206, 0.091920438733816695, 0.2020452320238027, 0.015783419431732466, 0.15438044913006899, 0.33626901196743414])
    #Use the estimates from the 1st run
    
    guess = np.log([0.49628073,  0.58854499 , 0.48857723,  2.10157473,  1.16726338])
    guess = np.append(guess,[0.547125332615])
    guess_w_blen = np.log([0.031386045599151487, 0.034027414206324061, 0.014031840798100775, 0.091769060449985645, 0.20503690355041698, 0.015704615057422228, 0.15445750717579823, 0.33510638411285476])
    #Use the estimates from the unfinished run
    guess_w_blen = np.append(guess_w_blen, guess)
    #r4 = test4.estimate(args, 7, guess_w_blen, est_blen = True,unrooted = False,force_tau = False)

    guess = np.log([0.49628073,  0.58854499 , 0.48857723,  2.10157473,  1.16726338])
    #guess = np.append(guess,[0.547125332615])
    guess_w_blen = np.log([0.031386045599151487, 0.034027414206324061, 0.014031840798100775, 0.091769060449985645, 0.20503690355041698, 0.015704615057422228, 0.15445750717579823, 0.33510638411285476])
    #Use the estimates from the unfinished run
    guess_w_blen = np.append(guess_w_blen, guess)
    
    r5 = test4.estimate(args, 7, guess_w_blen, est_blen = True,unrooted = False,force_tau = True)
    

##    #guess_w_blen[0:8]=[-1.0,0.0,1.0,-2.0,-3., -1., 0.0,-2.0]
##    guess_w_blen[0:2]=[-1.0,-0.3]
##    r2 = test2.estimate(args, sim_SubModel, guess_w_blen, True, True)
##
##    print 'blens are : ', np.exp(r1['x'][0:len(blen)])
##    print 'paras are : ', np.exp(r1['x'][len(blen):-1]),r1['x'][-1]
##    
##    print 'blens are : ', np.exp(r2['x'][0:len(blen)])
##    print 'paras are : ', np.exp(r2['x'][len(blen):-1]),r2['x'][-1]

    #cleanPAML('./mc.paml','./clnmc.paml',True)
        

        
