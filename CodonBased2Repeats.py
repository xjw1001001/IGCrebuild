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
import scipy.linalg
from numpy.testing import assert_allclose

import npmctree
from npmctree.sampling import sample_histories

from npmctree.sampling import sample_histories,sample_history
from npmctree.dynamic_lmap_lhood import get_iid_lhoods, get_lhood

from Bio import Phylo
from Bio import SeqIO
from cStringIO import StringIO
from Bio.Phylo import PhyloXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from copy import deepcopy



class Codon2RepeatsPhy:
    def __init__(self,numLeaf,blen,tree_newick,dataloc):
        self.nleaf = numLeaf
        self.nbranch = 2*numLeaf -2
        self.blen = blen
##        if not self.nbranch == len(blen):
##            print 'please make sure the vector blen contains all leaves'
        self.newicktree = tree_newick
        self.seqloc = dataloc
        self.nsites = 0
        self.treetopo, self.root, self.edge_to_blen, self.tree_phy = self.get_tree_info()
        self.models={'Jukes-Cantor':0,'K80':1, 'F81':2, 'HKY':3,'GTR':4, 'General Time-reversible':5,'Codon':6}
        self.duplication_node = 'N0' #This is for test version, need to change later
        self.SpecAfterDupli_node = 'N1' #This is for test version, need to change later 
        self.data = self.get_data()
        self.err = 1e-10
        
    def DataPre(self):
        cmdline = MuscleCommandline(input=self.seqloc, out="data_out.aln", clw=True)
        cmdline()
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

    def get_tree_info(self,align = False):
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
            lookup = dict((rec.id,rec.seq.tostring()) for rec in self.DataPre())
        else:
            fastaseqs = SeqIO.parse(open(self.seqloc,"rU"),'fasta')
            lookup = dict((rec.id,rec.seq.tostring()) for rec in fastaseqs)
            self.nsites = len(lookup[lookup.keys()[0]])
        
        for clade in tree_phy.get_terminals():
            key = clade.name            
            clade.sequences.append(lookup[key+'EDN'])
            if lookup.has_key(key+'ECP'):
                clade.sequences.append(lookup[key+'ECP'])
            else:
                print("Clade",key,"doesn't have ECP gene")
            
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name,v.name,d['weight']) for (u,v,d) in tree_nx.edges(data=True))
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

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

    def get_data(self):
        d = {}
        name = ['EDN','ECP']
        for clade in self.tree_phy.get_terminals():
            for i in range(len(clade.sequences)):
                d[clade.name,name[i]] = clade.sequences[i]
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
        bases = 'tcag'
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



        return Qnorm,dist
        

    def get_Q_and_distn(self,SubModel,para,Tao,NumDupli = 2):
        #By default, use Jukes-Cantor Model
        casenum = self.getcasenum(SubModel)
        Qnorm,dist = self.Q_norm(SubModel,para)
        expected_rate = np.dot(dist, -np.diag(Qnorm))
        Qnorm = Qnorm / expected_rate
        nt_pairs, pair_to_state = self.get_state_space(2)
        nstates = len(nt_pairs)
        n = len(nt_pairs)
        Q_un = np.zeros((n,n),dtype=float)
        if NumDupli == 1:
            if casenum<6:
                for i, (s0a, s1a) in enumerate(nt_pairs):
                    for j, (s0b, s1b) in enumerate(nt_pairs):
                        if i==j:
                            continue
                        if s0a == s1a and s0b == s1b:
                            Q_un[i,j] = Qnorm['ACGT'.index(s0a),'ACGT'.index(s0b)]
                Q_un = Q_un - np.diag(Q_un.sum(axis=1))
                distn = [0.0]*n#np.zeros((1,n),dtype=float)
                for i, (a,b) in enumerate(nt_pairs):
                    if a==b:
                        distn[i] = dist['ACGT'.index(a)]
                
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
            w,v = scipy.sparse.linalg.eigs(Q_un.T, k=1, which = 'SM')
            weights = v[:,0].real
            distn = weights / weights.sum()

            expected_rate = np.dot(distn, -np.diag(Q_un))
            Q = Q_un/expected_rate
            
#not normalized
            return Q_un,distn
        else:
            print('Please Check get_Q_and_distn function and make sure the case is considered')
            return Qnorm, distn
        
    def get_P(self,edge_to_blen_infer,SubModel,para,Tao):
        edge_to_P = {}
        for edge in self.treetopo.edges():
            blen = edge_to_blen_infer[edge]
            all_dupli_nodes = nx.descendants(self.treetopo,self.SpecAfterDupli_node)
            if edge[1] in all_dupli_nodes:
                Q,dist = self.get_Q_and_distn(SubModel,para,Tao,2)
            else:
                Q,distn = self.get_Q_and_distn(SubModel,para,Tao,1)
            P = scipy.sparse.linalg.expm(float(blen) * Q)
            edge_to_P[edge] = P
        root_distn = distn
        return edge_to_P, root_distn

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
        paralogs = ('EDN', 'ECP')
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
                nt0 = data[node, paralogs[0]][site]
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


    def estimate(self,args,SubModel,guess,est_blen = False):
##        self.treetopo, self.root, self.edge_to_blen, self.tree_phy = self.get_tree_info()
##        T, root, edge_to_blen, tree_phy = test.get_tree_info()
        nt_pairs, pair_to_state = self.get_state_space(2)
        nodes = set(self.treetopo)
        constraint_info = self.get_data_constraints(nodes, pair_to_state, self.nsites, self.data)
        constraints, npaired_yes, npaired_no = constraint_info
        
        #guess = [ 1.02444478 , 1.07097963 , 1.01565945 , 1.01566828 , 1.07064164 , 1.03295529,0.17774279,  0.33218723,  0.      ,    1.97231918 , 0.56235118]

        if est_blen:
            casenum = self.getcasenum(SubModel)
            one_bnd = (None,None)
            if casenum == 0:            
                bnds = [one_bnd]*(len(self.edge_to_blen)+2)
            elif casenum == 1:
                bnds = [one_bnd]*(len(self.edge_to_blen)+2)
            elif casenum == 2:
                bnds = [one_bnd]*len(self.edge_to_blen)
                bnds.extend([(0.0, 1.0),(0.0, 1.0),(0.0, 1.0),(0, None)])
            elif casenum == 3:
                bnds = [one_bnd]*len(self.edge_to_blen)
                bnds.extend([(0.0, 1.0),(0.0, 1.0),(0.0, 1.0),(0, None),(0, None)])

            ff = partial(self.objective_for_blen, SubModel, nt_pairs, constraints)
            
                     
        else:
            f = partial(self.objective, SubModel, self.treetopo, self.root, nt_pairs, constraints, self.edge_to_blen)
            casenum = self.getcasenum(SubModel)
            if casenum == 0:
                bnds = ((0, None),(0, None))                
            elif casenum == 1:
                bnds = ((0, None),(0, None))                
            elif casenum == 2:
                bnds = ((self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, None))                
            elif casenum == 3:
                bnds = [(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, 1.0-self.err),(self.err, None),(0, None)]                
            ff = lambda x:f(x[0:-1],x[-1])
##        f = partial(self.objective, SubModel, self.treetopo, self.root, self.edge_to_blen, nt_pairs, constraints,[1.0])
                
        print 'paralog site matches :', npaired_yes
        print 'paralog site mismatches:', npaired_no
        nsmall = 4
        pa = npaired_yes+0.0
        pb = npaired_no+0.0
        phi_hat = (pa+pb)/pb*3.0-4.0
        print 'preliminary estimate of phi:', guess,phi_hat
        print 'negative log likelihood of preliminary estimate:', ff(guess)        
        

        result = scipy.optimize.minimize(ff, x0=guess, method='L-BFGS-B', bounds=bnds)
        print result
        
    def objective_for_blen(self,SubModel,nt_pairs,constraints,x):
        casenum = self.getcasenum(SubModel)
        est_edge_to_blen = deepcopy(self.edge_to_blen)
        for i, k in enumerate(est_edge_to_blen):
            est_edge_to_blen[k] = np.exp(x[i])

        start_of_para = len(est_edge_to_blen)
        est_para = x[start_of_para:-1]
        est_Tao = x[-1]

        rst = self.objective(SubModel, self.treetopo, self.root, nt_pairs, constraints, est_edge_to_blen,est_para, est_Tao)

        print 'para:', x, rst

        #objective(self,SubModel, T, root,nt_pairs,constraints,edge_to_blen_infer,para,Tao)

        return rst


    def simulator(self,edge_to_blen_infer,SubModel,para,Tao,nsites,codon = False):
        edge_to_P, root_distn = self.get_P(edge_to_blen_infer,SubModel,para,Tao)
        node_to_data_lmap = {}
        if not codon:
            nt_pairs, pair_to_state = self.get_state_space(2)
            state_pairs = nt_pairs
        else:
            codon_pairs,pair_to_state = get_codon_state_space(2)
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

        
        
  

if __name__ == '__main__':
    numLeaf = 4
    blen = np.ones([2*numLeaf-2])*2
    tree_newick = './data/input_tree_test.newick'
    dataloc = './data/input_data.fasta'
    simdata = 'simdata.fasta'
    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc)
    
    sim_SubModel=3
    sim_para=[0.2,0.3,0.1,1.2]
    sim_Tao=1.0
    sim_nsites = 480
    guess = [0.4,0.3,0.6,2.0]
    guess.append(sim_Tao)
#    a = test.simulator(test.edge_to_blen,sim_SubModel,sim_para,sim_Tao,sim_nsites,False)
#    test.simtofasta(a,simdata,['EDN','ECP'])
    test2 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,simdata)
##    test.drawtree()

    parser = argparse.ArgumentParser()
##    parser.add_argument('--m', default=0, type=int,
##            help='submodel number, default 0 for Jukes-Cantor Model')
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()

##    test_main(args)

##    test.estimate(args,0)
#    test2.estimate(args,sim_SubModel,guess)
    guess_w_blen = [1.0]*len(test.edge_to_blen)
    guess_w_blen.extend(guess)
    test2.estimate(args, sim_SubModel, guess_w_blen, True)


        
