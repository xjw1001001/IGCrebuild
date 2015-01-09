# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonBased2Repeats class
# commit number: Oct 22nd, 2014 for old package
# cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5
# commit number: Dec 17th, 2014 for new package
# 33e393a973161e3a29149e82bfda23882b5826f3
from CodonGeneconFunc import *
import argparse
import ast

class CodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, Model = 'MG94', nnsites = None):
        self.newicktree = tree_newick  # newick tree file loc
        self.seqloc = alignment # multiple sequence alignment, now need to remove gap before-hand
        self.paralog = paralog  # paralog list
        self.nsites = nnsites  # number of sites in the alignment used for calculation
        self.Model = Model
        self.ll = 0.0


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
        self.get_data(Model)
         # self.name_to_seq (actrually name_to_codonlist)
         # self.observable_names
         # self.observable_nodes
         # self.observable_axes
         # self.iid_observations

        # Rate matrix related variable
        self.get_initial_x_process()
        self.update_by_x()
        self.unpack_x_process()
         # self.x_process  log values
         # self.x_rates    log values
         # self.x (x_process + x_rates)
         # self.x_clock (x_process + Lr)
         # self.pi         real values
         # self.kappa      real values
         # self.omega      real values
         # self.tau        real values
        self.processes = self.get_processes() #self.get_MG94Geneconv_and_MG94()
        self.GeneconvTransRed = None
        self.ExpectedGeneconv = None
        self.get_prior() 

        # guess related variable

    def get_processes(self):
        if self.Model == 'MG94':
            return self.get_MG94Geneconv_and_MG94()
        elif self.Model == 'HKY':
            return self.get_HKYGeneconv()

    def get_initial_x_process(self, kappa = 1.2, omega = 0.9, tau = 1.4, Force = None):
        count = np.array([0, 0, 0, 0], dtype = float) # count for A, C, G, T in all seq
        for name in self.name_to_seq:
            for i in range(4):
                count[i] += ''.join(self.name_to_seq[name]).count('ACGT'[i])
        count = count / count.sum()

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, omega, tau
            x_process = np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  kappa, omega, tau])
            x_process = np.log(x_process)
        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            x_process = np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  kappa, tau])
            x_process = np.log(x_process)            
        setattr(self, 'x_process', x_process)
        
        l = len(self.edge_to_blen) / 2 + 1               # number of leaves
        x_rates_clock = np.ones((l)) * 0.9
        x_clock = np.concatenate((x_process, np.log(x_rates_clock)))
        setattr(self, 'x_clock', x_clock)

        x_rates = np.array([ self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()])
        x_rates = np.log(x_rates)
        setattr(self, 'x_rates', x_rates)

        x = np.concatenate((x_process, x_rates))
        setattr(self, 'x', x)
        self.update_by_x_clock(Force = Force)

        if Force != None:
            if self.Model == 'MG94':
                self.update_by_x_clock(Force = Force)
            elif self.Model == 'HKY':
                self.update_by_x(Force = Force)

    def update_by_x(self, x = None, Force = None):
        k = len(self.edge_to_blen)
        if x != None:
            self.x = x
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        Force_process = None
        Force_rates = None
        if Force != None:
            Force_process = {i:Force[i] for i in Force.keys() if i < len(self.x) - k}
            Force_rates = {i:Force[i] for i in Force.keys() if not i < len(self.x) - k}
        self.unpack_x_process(Force = Force_process)
        self.unpack_x_rates(Force = Force_rates)

    def update_by_x_clock(self, x_clock = None, display = False, Force = None):
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        if x_clock != None:
            self.x_clock = x_clock
            
        self.x_process, Lr = self.x_clock[:-l], np.exp(self.x_clock[-l:])
        x_edge_clock = []

        for edge in self.edge_to_blen.keys():
            if edge[0] == 'N0':  # here I abondoned root node
                if str.isdigit(edge[1][1:]):  # (N0, N1) branch
                    x_edge_clock.append( Lr[0] * Lr[1] * (1 - Lr[2]) )
                else:
                    x_edge_clock.append( Lr[0] * (2 - Lr[1]) )

            else:
                tmp_k = int(edge[0][1:])
                if str.isdigit( edge[1][1:] ): # ( N_temp_k, N_temp_k+1 ) branch
                    x_edge_clock.append( reduce( mul, Lr[: (tmp_k + 2)], 1)  * (1 - Lr[tmp_k + 2]) )
                else:  # ( N_temp_k, leaf ) branch
                    x_edge_clock.append( reduce( mul, Lr[: (tmp_k + 2)], 1) )

        self.x_rates = np.log(np.array(x_edge_clock))
        self.x = np.concatenate([self.x_process, self.x_rates])
        self.update_by_x(Force = Force)  # TODO: Need to translate Force for blen constraints
        
        if display:
            print 'Lr rates: ', Lr
            print 'Process rates: ', np.exp(self.x_process)
            #print self.tau
        

    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        setattr(self, 'tree_phy', tree_phy)
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

        setattr(self, 'edge_to_blen', edge_to_blen)
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        Outgroup = list(set(leaves).difference(nx.descendants(T, self.SpecAfterDupli_node)))
        node_names = list(T)
        node_to_number = {n: i for i, n in enumerate(node_names)}
        setattr(self, 'node_to_num', node_to_number)
        tree_row = [self.node_to_num[na] for na, nb in edge_to_blen.keys()]
        tree_col = [self.node_to_num[nb] for na, nb in edge_to_blen.keys()]
        tree_process = [0 if e == ('N0', Outgroup[0]) else 1 for e in edge_to_blen]
        tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )
        setattr(self, 'tree', tree)

    def update_tree(self):
        tree_rate = [self.edge_to_blen[k] for k in self.edge_to_blen.keys()]
        self.tree['rate'] = tree_rate
    
    def unpack_x_process(self, Force = None):
        x_process = np.exp(self.x_process)

        if Force != None:
            for i in Force.keys():
                if i < len(self.x_process):
                    x_process[i] = Force[i]

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, tau, omega
            assert(len(self.x_process) == 6)
            
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.omega = x_process[4]
            self.tau = x_process[5]
        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            assert(len(self.x_process) == 5)
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.tau = x_process[4]




        # Now update the prior distribution
        self.get_prior()

        # Now update processes (Rate matrices)
        self.get_processes()
        

    def unpack_x_rates(self, Force = None):  # TODO: Change it to fit general tree structure rather than cherry tree
        x_rates = np.exp(self.x_rates)
        if Force != None:
            for i in Force.keys():
                if i < len(self.x_rates):
                    x_rates[i] = Force[i]
        assert(len(x_rates) == len(self.edge_to_blen))

        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        # Now update blen with fixed order:
        # Always start from the root and internal-tip branch first
        for i in range(len(internal_branch)):
            self.edge_to_blen[internal_branch[i]] = x_rates[2 * i]
            self.edge_to_blen[leaf_branch[i]] = x_rates[2 * i + 1]
        for j in range(len(leaf_branch[i + 1:])):
            self.edge_to_blen[leaf_branch[i + 1 + j]] = x_rates[ - len(leaf_branch[i + 1:]) + j]
##        self.edge_to_blen = {}
##        for it in range(len(self.edge_to_blen)):
##            edge = self.edge_to_blen.keys()[it]
##            self.edge_to_blen[edge] = x_rates[it]
        self.update_tree()
         

    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(self.name_to_seq[name]) / 3 )]
            self.name_to_seq[name] = tmp_seq
            
        
    def get_data(self, Model):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        setattr(self, 'name_to_seq', {name:str(seq_dict[name].seq) for name in seq_dict.keys()})

        if Model == 'MG94':
            # Convert from nucleotide sequences to codon sequences.
            self.nts_to_codons()
            obs_to_state = self.codon_to_state
        else:
            obs_to_state = self.nt_to_state
        
        # Throttle the number of sites if requested.
        if self.nsites is None:
            self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])
        else:
            for name in self.name_to_seq:
                self.name_to_seq[name] = self.name_to_seq[name][: self.nsites]

        print 'number of sites to be analyzed: ', self.nsites

        observable_names = self.name_to_seq.keys()
        setattr(self, 'observable_names', observable_names)
        paralog_len = [len(a) for a in self.paralog]
        assert(paralog_len[1:] == paralog_len[:-1])  # check if all paralog names have same length
        suffix_len = len(self.paralog[0])
        #observable_suffixes = [name[-suffix_len:] for name in observable_names]
        observable_suffixes = self.paralog
        suffix_to_axis = {n:i for (i, n) in enumerate(list(set(observable_suffixes))) }
        observable_nodes = [self.node_to_num[n[:-suffix_len]] for n in observable_names]
        observable_axes = [suffix_to_axis[s[-suffix_len:]] for s in observable_names]
        setattr(self, 'observable_axes', observable_axes)
        setattr(self, 'observable_nodes', observable_nodes)


        iid_observations = []
        for site in range(self.nsites):
            observations = []
            for name in self.observable_names:
                observation = obs_to_state[self.name_to_seq[name][site]]
                observations.append(observation)
            iid_observations.append(observations)
        setattr(self, 'iid_observations', iid_observations)

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
        if self.Model == 'MG94':
            prior_feasible_states = [(self.codon_to_state[codon], self.codon_to_state[codon]) for codon in self.codon_nonstop]
            distn = [ reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
            distn = np.array(distn) / sum(distn)
        elif self.Model == 'HKY':
            prior_feasible_states = [(self.nt_to_state[nt], self.nt_to_state[nt]) for nt in 'ACGT']
            distn = [ self.pi['ACGT'.index(nt)] for nt in 'ACGT' ]
            distn = np.array(distn) / sum(distn)
        setattr(self, 'prior_feasible_states', prior_feasible_states)
        setattr(self, 'prior_distribution', distn)

    def vec_get_MG94Geneconv(self):
        pair_list = [a + b for a, b in product(self.codon_nonstop, repeat = 2)]
        pair_from_matrix = np.repeat(pair_list, len(pair_list)).reshape((len(pair_list), len(pair_list)))
        process_geneconv = vec_get_MG94GeneconvRate(pair_from = pair_from_matrix, pair_to = pair_list, pi = self.pi,
                                   kappa = self.kappa, omega = self.omega, codon_table = self.codon_table,
                                   tau = self.tau, codon_to_state = self.codon_to_state)
    
    def get_HKYGeneconv(self, x = None):
        #print 'tau = ', self.tau
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        #self.update_by_x(x)

        for i, pair_from in enumerate(product('ACGT', repeat = 2)):
            na, nb = pair_from
            sa = self.nt_to_state[na]
            sb = self.nt_to_state[nb]
            for j, pair_to in enumerate(product('ACGT', repeat = 2)):
                nc, nd = pair_to
                sc = self.nt_to_state[nc]
                sd = self.nt_to_state[nd]
                if i == j:
                    continue
                GeneconvRate = get_HKYGeneconvRate(pair_from, pair_to, self.pi, self.kappa, self.tau)
                if GeneconvRate != 0.0:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(0.0)
                if na == nb and nc == nd:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(get_HKYBasicRate(na, nc, self.pi, self.kappa))

        process_geneconv = dict(
            row = row,
            col = col,
            rate = rate_geneconv
            )
        process_basic = dict(
            row = row,
            col = col,
            rate = rate_basic
            )
        return [process_basic, process_geneconv]
                    
            


    def get_MG94Geneconv_and_MG94(self, x = None):
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        #self.update_by_x(x)

        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
            #rate_sum_geneconv = 0
            #rate_sum_basic = 0
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
                    Qbasic = get_MG94BasicRate(cb, cc, self.pi, self.kappa, self.omega, self.codon_table)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qbasic)
                        #rate_sum_geneconv += Qbasic
                        rate_basic.append(0.0)
                        #rate_sum_basic += 0.0

                    # (ca, cb) to (cc, cb)
                    Qbasic = get_MG94BasicRate(ca, cc, self.pi, self.kappa, self.omega, self.codon_table)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate_geneconv.append(Qbasic)
                        #rate_sum_geneconv += Qbasic
                        rate_basic.append(0.0)
                        #rate_sum_basic += 0.0

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qbasic = get_MG94BasicRate(cb, ca, self.pi, self.kappa, self.omega, self.codon_table)
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                rate_geneconv.append(Qbasic + Tgeneconv)
                #rate_sum_geneconv += Qbasic + Tgeneconv
                rate_basic.append(0.0)
                #rate_sum_basic += 0.0
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qbasic = get_MG94BasicRate(ca, cb, self.pi, self.kappa, self.omega, self.codon_table)
                rate_geneconv.append(Qbasic + Tgeneconv)
                #rate_sum_geneconv += Qbasic + Tgeneconv
                rate_basic.append(0.0)
                #rate_sum_basic += 0.0

                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                #row.append((sa, sb))
                #col.append((sa, sb))
                #rate_geneconv.append(-rate_sum_geneconv)
                #rate_basic.append(-rate_sum_basic)
            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qbasic = get_MG94BasicRate(ca, cc, self.pi, self.kappa, self.omega, self.codon_table)
                    if Qbasic != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qbasic)
                        rate_basic.append(0.0)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate_geneconv.append(Qbasic)
                        rate_basic.append(0.0)
                        #rate_sum_geneconv += 2 * Qbasic
                        #rate_sum_basic += 2 * Qbasic

                    # (ca, ca) to (cc, cc)
                        row.append((sa, sb))
                        col.append((sc, sc))
                        rate_geneconv.append(0.0)
                        rate_basic.append(Qbasic)
                        
                # Finally add the diagonal
                # (ca, cb) to (ca, cb)
                #row.append((sa, sb))
                #col.append((sa, sb))
                #rate_geneconv.append(-rate_sum_geneconv)
                #rate_basic.append(-rate_sum_basic)
                
        process_geneconv = dict(
            row = row,
            col = col,
            rate = rate_geneconv
            )
        process_basic = dict(
            row = row,
            col = col,
            rate = rate_basic
            )
        return [process_basic, process_geneconv]

    def _loglikelihood(self, edge_derivative = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
##        k = len(self.edge_to_blen)
##        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        if self.Model == 'MG94':
##            self.processes = self.get_MG94Geneconv_and_MG94()
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
##            self.processes = self.get_HKYGeneconv()
            state_space_shape = [4, 4]

        # prepare some extra parameters for the json interface
        if edge_derivative:
            requested_derivatives = list(range(k))
        else:
            requested_derivatives = []
            
        site_weights = np.ones(self.nsites)

        # prepare the input for the json interface
        data = dict(
            site_weights = site_weights,
            requested_derivatives = requested_derivatives,
            node_count = len(self.edge_to_blen) + 1,
            state_space_shape = state_space_shape,
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

        status = j_ll['status']
        feasibility = j_ll['feasibility']

        if status != 'success' or not feasibility:
            print 'results:'
            print j_ll
            print
            raise Exception('Encountered some problem in the calculation of log likelihood and its derivatives')

        ll, edge_derivs = j_ll['log_likelihood'], j_ll['edge_derivatives']

        return ll, edge_derivs

    def _loglikelihood2(self, edge_derivative = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        scene = self.get_scene()
        
        log_likelihood_request = {'property':'snnlogl'}
        derivatives_request = {'property':'sdnderi'}
        if edge_derivative:
            requests = [log_likelihood_request, derivatives_request]
        else:
            requests = [log_likelihood_request]
        j_in = {
            'scene' : scene,
            'requests' : requests
            }
        j_out = jsonctmctree.interface.process_json_in(j_in)

        status = j_out['status']
    
        ll = j_out['responses'][0]
        if edge_derivative:
            edge_derivs = j_out['responses'][1]
        else:
            edge_derivs = []

        return ll, edge_derivs
    
    def get_scene(self):
        #k = len(self.edge_to_blen)
        #self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        if self.Model == 'MG94':
            #self.processes = self.get_MG94Geneconv_and_MG94()
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
            #self.processes = self.get_HKYGeneconv()
            state_space_shape = [4, 4]
        process_definitions = [{'row_states':i['row'], 'column_states':i['col'], 'transition_rates':i['rate']} for i in self.processes]
        scene = dict(
            node_count = len(self.edge_to_blen) + 1,
            process_count = len(self.processes),
            state_space_shape = state_space_shape,
            tree = {
                'row_nodes' : self.tree['row'],
                'column_nodes' : self.tree['col'],
                'edge_rate_scaling_factors' : self.tree['rate'],
                'edge_processes' : self.tree['process']
                },
            root_prior = {'states':self.prior_feasible_states,
                          'probabilities': self.prior_distribution},
            process_definitions = process_definitions,
            observed_data = {
                'nodes':self.observable_nodes,
                'variables':self.observable_axes,
                'iid_observations':self.iid_observations
                }            
            )
        return scene
    
    def get_geneconvTransRed(self, get_rate = False):
        row_states = []
        column_states = []
        if self.Model == 'MG94':
            for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
                ca, cb = pair
                sa = self.codon_to_state[ca]
                sb = self.codon_to_state[cb]
                if ca == cb:
                    continue
                
                # (ca, cb) to (ca, ca)
                row_states.append((sa, sb))
                column_states.append((sa, sa))
                # (ca, cb) to (cb, cb)
                row_states.append((sa, sb))
                column_states.append((sb, sb))
            
        elif self.Model == 'HKY':
            for i, pair in enumerate(product('ACGT', repeat = 2)):
                na, nb = pair
                sa = self.nt_to_state[na]
                sb = self.nt_to_state[nb]
                if na == nb:
                    continue

                # (na, nb) to (na, na)
                row_states.append((sa, sb))
                column_states.append((sa, sa))

                # (na, nb) to (nb, nb)
                row_states.append((sa, sb))
                column_states.append((sb, sb))

        return {'row_states' : row_states, 'column_states' : column_states, 'weights' : np.ones(len(row_states))}
                
                


    def _ExpectedNumGeneconv(self, package = 'new', display = False):
        if self.GeneconvTransRed is None:
            self.GeneconvTransRed = self.get_geneconvTransRed()

        if package == 'new':
            scene = self.get_scene()
            requests = [{'property' : 'SDNTRAN', 'transition_reduction' : self.GeneconvTransRed}]
            j_in = {
                'scene' : scene,
                'requests' : requests
                }        
            j_out = jsonctmctree.interface.process_json_in(j_in)

            status = j_out['status']
            ExpectedGeneconv = {self.edge_to_blen.keys()[i] : j_out['responses'][0][i] for i in range(len(self.edge_to_blen))}
            return ExpectedGeneconv
        else:
            print 'Need to implement this for old package'

    def get_ExpectedNumGeneconv(self):
        self.ExpectedGeneconv = self._ExpectedNumGeneconv()
            
    def loglikelihood_and_gradient(self, package = 'new', display = False, Force = None):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        self.update_by_x(Force = Force)
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array
        if package == 'new':
            fn = self._loglikelihood2
        else:
            fn = self._loglikelihood

        ll, edge_derivs = fn(edge_derivative = True)
        
        m = len(self.x) - len(self.edge_to_blen)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
            if Force != None:
                if i in Force.keys():
                    other_derivs.append(0.0)
                    continue
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta
            self.update_by_x(x_plus_delta, Force = Force)
            ll_delta, _ = fn(edge_derivative = False)
            d_estimate = (ll_delta - ll) / delta           
            other_derivs.append(d_estimate)
            # restore self.x
            self.update_by_x(x, Force = Force)
        other_derivs = np.array(other_derivs)
        if display:
            print 'log likelihood = ', ll
            print 'Edge derivatives = ', edge_derivs
            print 'other derivatives:', other_derivs
            print 'Current x array = ', self.x

        self.ll = ll
        f = -ll
        g = -np.concatenate((other_derivs, edge_derivs))
        return f, g
        #ll, edge_derivs =

    def Clock_wrap(self, display, Force, x_clock):
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        self.update_by_x_clock(x_clock, display, Force = Force)
        self.x_process, Lr = x_clock[:-l], np.exp(x_clock[-l:])

        f, g = self.loglikelihood_and_gradient(Force = Force)
        
        # Now need to calculate the derivates
        other_derives, edge_derives = g[:-nEdge], g[-nEdge:]
        edge_to_derives = {self.edge_to_blen.keys()[i] : edge_derives[i] for i in range(len(self.edge_to_blen.keys()))}

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        Lr_derives = []  # used to store derivatives for the clock parameters L, r0, r1, ...
        Lr_derives.append(sum(edge_derives))  # dLL/dL = sum(all derives)
        Lr_derives.append(edge_to_derives[out_group_branch] * 2 / (Lr[1] - 2)
                          + sum(edge_derives))

        for i in range(2, len(Lr)):  # r(i-1)
            Lr_derives.append( edge_to_derives[('N' + str(i - 2), 'N' + str(i - 1))] * Lr[i] / (Lr[i] - 1) # 
                               + sum([edge_to_derives[internal_branch[j]] for j in range(i - 1, len(internal_branch))])  # only sum over nodes decendent from node i-1
                               + sum([edge_to_derives[leaf_branch[j]] for j in range(i - 1, len(leaf_branch))]))  # only sum over nodes decendent from node i-1

        #TODO: Need to change the two sums if using general tree

        g_clock = np.concatenate( (np.array(other_derives), np.array(Lr_derives)))

        if display:
            print 'log likelihood = ', f
            print 'Lr derivatives = ', Lr_derives
            print 'other derivatives = ', other_derives
            print 'Current x_clock array = ', self.x_clock

        return f, g_clock
        
    def objective_and_gradient(self, display, Force, x):
        self.update_by_x(x, Force = Force)
        f, g = self.loglikelihood_and_gradient(display = display, Force = Force)
        return f, g

    def get_mle(self, clock = False, display = True, Force = None):
        bnds = [(None, -0.05)] * 3
        if not clock:
            self.update_by_x(Force = Force)
            f = partial(self.objective_and_gradient, display, Force)
            guess_x = self.x            
            bnds.extend([(None, None)] * (len(self.x) - 3))
            
        else:
            self.update_by_x_clock(Force = Force)  # TODO: change force for blen in x_clock
            f = partial(self.Clock_wrap, display, Force)
            guess_x = self.x_clock
            bnds.extend([(None, None)] * (len(self.x_clock) - 2 - (len(self.edge_to_blen) / 2 + 1)))
            bnds.extend([(None, 0.0)] * (len(self.edge_to_blen) / 2))        

        result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        print (result)
        return result

    def save_to_file(self, clock = False, file_name = None, path = './', Force = None):
        if file_name == None:
            file_name = self.Model + '_' + '_'.join(self.paralog)
            if clock:
                file_name += '_clock.p'
            else:
                file_name += '_nonclock.p'
        save_file = path + file_name
        self.ll = self._loglikelihood()[0]
        print 'x = ', self.x, 'x_clock = ', self.x_clock
        save_info = dict(
            Model = self.Model,
            x     = self.x.tolist(),
            x_clock = self.x_clock.tolist(),
            edge_to_blen = self.edge_to_blen,
            ExpectedGeneconv = self.ExpectedGeneconv,
            ll    = self.ll,
            newicktree = self.newicktree,
            alignment_file = self.seqloc,
            paralog = self.paralog,
            clock = clock,
            Force = Force
            )
        pickle.dump(save_info, open(save_file, 'wb+'))  # use pickle to save the class which can be easily reconstructed by pickle.load()
        
def main(args):
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../PairsAlignemt/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94')
    path = './NewPackageNewRun/'

    print 'Now calculate MLE for pair', paralog
    x = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574, 1.19617746,
                   -3.90910406, -2.01786889,  # N0
                   -2.92359461, -2.83555499,  # N1
                   -4.30005794, -3.29056063,  # N2
                   -4.59465356, -3.49891453,  # N3
                   -6.28014701, -3.10701318,  # N4
                   -3.87159909, -3.97438916]) # N5
    
    if hasattr(args, 'Force'):
        Force = args.Force
        if Force != None:
            path += 'Force_'
    else:
        Force = None
    test.update_by_x(x, Force = Force)
    
    result = test.get_mle(clock = False, display = True, Force = Force)
    test.get_ExpectedNumGeneconv()
    test.save_to_file(clock = False, path = path, Force = Force)

    x_clock = np.array([-0.64430238, -0.48627304, -0.93852343, 1.76793238, -2.45061608, 1.23836258,
                  -2.11088631, -0.15313429, -0.23202426, -0.63091448, -0.25441366, -0.20629961,-0.23301088])
    test.update_by_x_clock(x_clock, Force = Force)
    result = test.get_mle(clock = True, display = True, Force = Force)
    test.get_ExpectedNumGeneconv()
    test.save_to_file(clock = True, path = path, Force = Force)
    

if __name__ == '__main__':##    parser = argparse.ArgumentParser()
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--Force', type = ast.literal_eval, help = 'Parameter constraints')
##    
##    main(parser.parse_args())

    paralog1 = 'YML026C'
    paralog2 = 'YDR450W'
    paralog1 = 'YDR502C'
    paralog2 = 'YLR180W'
    Force    = {5:0.0}
   # Force    = {4:0.0}
    #Force = None

    paralog = [paralog1, paralog2]
    alignment_file = '../PairsAlignemt/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94')#, nnsites = 1)
    
##    x = [ -0.71504098,  -0.51228865,  -0.89531866,   1.58369229,
##        -2.83087967,   1.19617746,  -2.82848596,  -3.37970062,
##        -4.27815837,  -4.98465689,  -2.17551306,  -3.22866221,
##       -15.46376865,  -2.73078175,  -4.01541511,  -1.91221709,
##       -18.25127742,  -2.62498074]
    
    xx = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574,
        1.19617746, -3.90910406, -3.87159909, -4.30005794, -6.28014701,
       -2.83555499, -3.29056063, -2.92359461, -3.49891453, -3.97438916,
       -3.10701318, -4.59465356, -2.01786889])

    #xx = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574,
    xx = np.array([-0.60829003, -0.51797889, -0.8170342,  1.76977596, -2.44028574, 1.19617746,
    #xx = np.array([-0.60829003, -0.51797889, -0.8170342,  1.76977596, -2.44028574,
                   -3.90910406, -2.01786889,  # N0
                   -2.92359461, -2.83555499,  # N1
                   -4.30005794, -3.29056063,  # N2
                   -4.59465356, -3.49891453,  # N3
                   -6.28014701, -3.10701318,  # N4
                   -3.87159909, -3.97438916]) # N5

    x = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574, 1.19617746,
    #x = np.array([-0.60829003, -0.51797889, -0.8170342,  1.76977596, -2.44028574,
                   -3.90910406, -2.01786889,  # N0
                   -2.92359461, -2.83555499,  # N1
                   -4.30005794, -3.29056063,  # N2
                   -4.59465356, -3.49891453,  # N3
                   -6.28014701, -3.10701318,  # N4
                   -3.87159909, -3.97438916]) # N5

    
    #test.update_by_x(x = xx, Force = Force)
    test.update_by_x(x = x)
    
    
    print 'tau = ', test.tau
    #print test._loglikelihood2()
    
    branch = test.edge_to_blen.keys()
    branch.sort(key = lambda node: int(node[0][1:]))
    for i in branch:
        print i, test.edge_to_blen[i]
##
##    tree_branch = [(test.tree['row'][i], test.tree['col'][i], test.tree['rate'][i]) for i in range(len(test.tree['col']))]
##    tree_branch.sort(key = lambda node: int(node[0]))
##    for i in tree_branch:
##        print i
##
##    print test.node_to_num

##    for i in range(20, 30):
##        print test.processes[0]['row'][i], test.processes[0]['col'][i], test.processes[0]['rate'][i], test.processes[1]['row'][i], test.processes[1]['col'][i], test.processes[1]['rate'][i]

        
