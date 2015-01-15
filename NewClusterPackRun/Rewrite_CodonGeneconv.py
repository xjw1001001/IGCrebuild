# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-write of my previous CodonGeneconv class
# commit number: Oct 22nd, 2014 for old package
# cb1ba60ee2b57d6703cd9a3987000c2fd4dd68a5
# commit number: Dec 17th, 2014 for new package
# 33e393a973161e3a29149e82bfda23882b5826f3

from CodonGeneconFunc import *
import argparse
import ast

class ReCodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, Model = 'MG94', nnsites = None, clock = False, Force = None):
        self.newicktree  = tree_newick  # newick tree file loc
        self.seqloc      = alignment    # multiple sequence alignment, now need to remove gap before-hand
        self.paralog     = paralog      # parlaog list
        self.nsites      = nnsites      # number of sites in the alignment used for calculation
        self.Model       = Model
        self.ll          = 0.0          # Store current log-likelihood
        self.Force       = Force        # parameter constraints only works on self.x not on x_clock which should be translated into self.x first
        self.clock       = clock        # molecular clock control

        # Tree topology related variable
        self.tree         = None        # store the tree dictionary used for json likelihood package parsing
        self.edge_to_blen = None        # dictionary store the unpacked tree branch length information {(node_from, node_to):blen}
        self.edge_list    = None        # kept all edges in the same order with x_rates
        self.node_to_num  = None        # dictionary used for translating tree info from self.edge_to_blen to self.tree
        self.num_to_node  = None        # dictionary used for translating tree info from self.tree to self.edge_to_blen

        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}

        # Tip data related variable
        self.name_to_seq      = None    # dictionary store sequences
        self.observable_names = None    # list of extent species + paralog name ( = self.name_to_seq.keys())
        self.observable_nodes = None    # list of extent species numbers (node_to_num)
        self.observable_axes  = None    # list of paralog numbers  
        self.iid_observations = None    # list of multivariate states


        # Rate matrix related variable
        self.x_process      = None      # log values of process parameters
        self.x_rates        = None      # log values of blen (general for both clock and nonclock cases)
        self.x              = None      # x_process + x_rates
        self.x_Lr           = None      # log values of clock blen parameters
        self.x_clock        = None      # x_process + Lr
        self.pi             = None      # real values
        self.kappa          = 1.2       # real values
        self.omega          = 0.9       # real values
        self.tau            = 1.4       # real values

        self.processes      = None      # list of basic and geneconv rate matrices. Each matrix is a dictionary used for json parsing

        # Prior distribution on the root
        self.prior_feasible_states  = None
        self.prior_distribution     = None

        # Expected_Geneconv Events
        self.GeneconvTransRed = None    # dictionary of Geneconv transition matrix used for json parsing
        self.ExpectedGeneconv = None    # dictionary storing expected number of geneconv events on each branch

        # Initialize all parameters
        self.initialize_parameters()
        
        
    def initialize_parameters(self):
        self.get_tree()
        self.get_data()
        self.get_initial_x_process()

    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        self.edge_to_blen = edge_to_blen

        # Now assign node_to_num
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        internal_nodes = set(list(T)).difference(leaves)
        node_names = list(internal_nodes) + list(leaves)
        self.node_to_num = {n:i for i, n in enumerate(node_names)}
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}

        # Prepare for generating self.tree so that it has same order as the self.x_process
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        edge_list = []
        for i in range(len(internal_branch)):
            edge_list.append(internal_branch[i])
            edge_list.append(leaf_branch[i])
        for j in range(len(leaf_branch[i + 1:])):
            edge_list.append(leaf_branch[i + 1 + j])
            
        # Now setup self.tree dictionary
        tree_row = [self.node_to_num[na] for na, nb in edge_list]
        tree_col = [self.node_to_num[nb] for na, nb in edge_list]
        tree_process = [0 if e[0] == 'N0' and e[1] != 'N1' else 1 for e in edge_list]
        self.edge_list = edge_list

        self.tree = dict(
            row = tree_row,
            col = tree_col,
            process = tree_process,
            rate = np.ones(len(tree_row))
            )

    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(self.name_to_seq[name]) / 3 )]
            self.name_to_seq[name] = tmp_seq
       
    def get_data(self):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        self.name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
        
        if self.Model == 'MG94':
            # Convert from nucleotide sequences to codon sequences.
            self.nts_to_codons()
            obs_to_state = self.codon_to_state
        else:
            obs_to_state = self.nt_to_state        

        # change the number of sites for calculation if requested
        if self.nsites is None:
            self.nsites = len(self.name_to_seq[self.name_to_seq.keys()[0]])
        else:
            for name in self.name_to_seq:
                self.name_to_seq[name] = self.name_to_seq[name][: self.nsites]
        print 'number of sites to be analyzed: ', self.nsites

        # assign observable parameters
        self.observable_names = self.name_to_seq.keys()
        paralog_len = [len(a) for a in self.paralog]
        assert(paralog_len[1:] == paralog_len[:-1])  # check if all paralog names have same length
        suffix_len = len(self.paralog[0])
        suffix_to_axis = {n:i for (i, n) in enumerate(list(set(self.paralog))) }
        self.observable_nodes = [self.node_to_num[n[:-suffix_len]] for n in self.observable_names]
        self.observable_axes = [suffix_to_axis[s[-suffix_len:]] for s in self.observable_names]
        
        # Now convert alignment into state list
        iid_observations = []
        for site in range(self.nsites):
            observations = []
            for name in self.observable_names:
                observation = obs_to_state[self.name_to_seq[name][site]]
                observations.append(observation)
            iid_observations.append(observations)
        self.iid_observations = iid_observations

    def get_initial_x_process(self):
        
        count = np.array([0, 0, 0, 0], dtype = float) # count for A, C, G, T in all seq
        for name in self.name_to_seq:
            for i in range(4):
                count[i] += ''.join(self.name_to_seq[name]).count('ACGT'[i])
        count = count / count.sum()

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, omega, tau
            self.x_process = np.log(np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  self.kappa, self.omega, self.tau]))
        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            self.omega = 1.0
            self.x_process = np.log(np.array([count[0] + count[2], count[0] / (count[0] + count[2]), count[1] / (count[1] + count[3]),
                                  self.kappa, self.tau]))

        self.x_rates = np.log(np.array([ self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()]))
        self.x = np.concatenate((self.x_process, self.x_rates))


        if self.clock:   # set-up x_clock if it's a clock model
            l = len(self.edge_to_blen) / 2 + 1               # number of leaves
            self.x_Lr = np.log(np.ones((l)) * 0.9)
            self.x_clock = np.concatenate((self.x_process, self.x_Lr))
            self.unpack_x_clock()

        self.update_by_x()
        
    def update_by_x_clock(self, x_clock = None):
        if not x_clock == None:
            self.x_clock = x_clock
        self.unpack_x_clock()
        self.update_by_x()
        
    def unpack_x_clock(self):
        assert(self.clock)
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.
        self.x_process, self.x_Lr = self.x_clock[:-l], np.exp(self.x_clock[-l:])

        # Now update self.x by using self.x_clock
        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        # Now update blen with fixed order:
        # Always start from the root and internal-tip branch first
        for i in range(len(internal_branch)):
            edge = internal_branch[i]
            self.x_rates[2 * i] = np.log(self.blen_from_clock(edge))
            edge = leaf_branch[i]
            self.x_rates[2 * i + 1] = np.log(self.blen_from_clock(edge))
        for j in range(len(leaf_branch[i + 1:])):
            edge = leaf_branch[i + 1 + j]
            self.x_rates[ - len(leaf_branch[i + 1:]) + j] = np.log(self.blen_from_clock(edge))
        # update self.x so that all parameters can be updated by update_by_x
        self.x = np.concatenate((self.x_process, self.x_rates))
        
    def blen_from_clock(self, edge):
        assert(edge in self.edge_to_blen.keys())
        if edge[0] == 'N0':
            if str.isdigit(edge[1][1:]):  # (N0, N1) branch
                return self.x_Lr[0] * self.x_Lr[1] * (1 - self.x_Lr[2])
            else:
                return self.x_Lr[0] * (2 - self.x_Lr[1])

        else:
            tmp_k = int(edge[0][1:])
            if str.isdigit( edge[1][1:] ): # ( N_temp_k, N_temp_k+1 ) branch
                return reduce( mul, self.x_Lr[: (tmp_k + 2)], 1)  * (1 - self.x_Lr[tmp_k + 2])
            else:  # ( N_temp_k, leaf ) branch
                return reduce( mul, self.x_Lr[: (tmp_k + 2)], 1)
        
        
    def update_by_x(self, x = None):
        k = len(self.edge_to_blen)
        if x != None:
            self.x = x
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        Force_process = None
        Force_rates = None
        if self.Force != None:
            Force_process = {i:self.Force[i] for i in self.Force.keys() if i < len(self.x) - k}
            Force_rates = {(i-k):self.Force[i] for i in self.Force.keys() if not i < len(self.x) - k}
        self.unpack_x_process(Force_process = Force_process)
        self.unpack_x_rates(Force_rates = Force_rates)

    def unpack_x_process(self, Force_process = None):
        x_process = np.exp(self.x_process)

        if Force_process != None:
            for i in Force_process.keys():
                x_process[i] = Force_process[i]

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

    def get_prior(self):
        if self.Model == 'MG94':
            self.prior_feasible_states = [(self.codon_to_state[codon], self.codon_to_state[codon]) for codon in self.codon_nonstop]
            distn = [ reduce(mul, [self.pi['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
            distn = np.array(distn) / sum(distn)
        elif self.Model == 'HKY':
            self.prior_feasible_states = [(self.nt_to_state[nt], self.nt_to_state[nt]) for nt in 'ACGT']
            distn = [ self.pi['ACGT'.index(nt)] for nt in 'ACGT' ]
            distn = np.array(distn) / sum(distn)
        self.prior_distribution = distn

    def get_processes(self):
        if self.Model == 'MG94':
            self.processes = self.get_MG94Geneconv_and_MG94()
        elif self.Model == 'HKY':
            self.processes = self.get_HKYGeneconv()

    def get_MG94Geneconv_and_MG94(self):
        Qbasic = self.get_MG94Basic()
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        for i, pair in enumerate(product(self.codon_nonstop, repeat = 2)):
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
                    Qb = Qbasic[sb, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                    # (ca, cb) to (cc, cb)
                    Qb = Qbasic[sa, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qb = Qbasic[sb, sa]
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv = self.tau * self.omega
                else:
                    Tgeneconv = self.tau
                rate_geneconv.append(Qb + Tgeneconv)
                rate_basic.append(0.0)
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qb = Qbasic[sa, sb]
                rate_geneconv.append(Qb + Tgeneconv)
                rate_basic.append(0.0)

            else:
                for cc in self.codon_nonstop:
                    if cc == ca:
                        continue
                    sc = self.codon_to_state[cc]

                    # (ca, ca) to (ca,  cc)
                    Qb = Qbasic[sa, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sa, sc))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)

                    # (ca, ca) to (cc, cc)
                        row.append((sa, sb))
                        col.append((sc, sc))
                        rate_geneconv.append(0.0)
                        rate_basic.append(Qb)
                
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

    def get_MG94Basic(self):
        Qbasic = np.zeros((61, 61), dtype = float)
        for ca in self.codon_nonstop:
            for cb in self.codon_nonstop:
                if ca == cb:
                    continue
                Qbasic[self.codon_to_state[ca], self.codon_to_state[cb]] = get_MG94BasicRate(ca, cb, self.pi, self.kappa, self.omega, self.codon_table)
        expected_rate = np.dot(self.prior_distribution, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic

    def get_HKYBasic(self):
        Qbasic = np.array([
            [0, 1.0, self.kappa, 1.0],
            [1.0, 0, 1.0, self.kappa],
            [self.kappa, 1.0, 0, 1.0],
            [1.0, self.kappa, 1.0, 0],
            ]) * np.array(self.pi)
        expected_rate = np.dot(self.prior_distribution, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic
        
    
    def get_HKYGeneconv(self):
        #print 'tau = ', self.tau
        Qbasic = self.get_HKYBasic()
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

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
                GeneconvRate = get_HKYGeneconvRate(pair_from, pair_to, Qbasic, self.tau)
                if GeneconvRate != 0.0:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(0.0)
                if na == nb and nc == nd:
                    row.append((sa, sb))
                    col.append((sc, sd))
                    rate_geneconv.append(GeneconvRate)
                    rate_basic.append(Qbasic['ACGT'.index(na), 'ACGT'.index(nc)])

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
    
    def unpack_x_rates(self, Force_rates = None):  # TODO: Change it to fit general tree structure rather than cherry tree
        x_rates = np.exp(self.x_rates)
        if Force_rates != None:
            for i in Force_rates.keys():
                x_rates[i] = Force_rates[i]
        assert(len(x_rates) == len(self.edge_to_blen))

        for edge_it in range(len(self.edge_list)):
            self.edge_to_blen[self.edge_list[edge_it]] = x_rates[edge_it] 

        self.update_tree()

    def update_tree(self):
        for i in range(len(self.tree['rate'])):
            node1 = self.num_to_node[self.tree['row'][i]]
            node2 = self.num_to_node[self.tree['col'][i]]
            self.tree['rate'][i] = self.edge_to_blen[(node1, node2)]

    def _loglikelihood(self, edge_derivative = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        if self.Model == 'MG94':
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
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
        if self.Model == 'MG94':
            state_space_shape = [61, 61]
        elif self.Model == 'HKY':
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

    def loglikelihood_and_gradient(self, package = 'new', display = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        self.update_by_x()
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
            if self.Force != None:
                if i in self.Force.keys():  # check here
                    other_derivs.append(0.0)
                    continue
            x_plus_delta = np.array(self.x)
            x_plus_delta[i] += delta
            self.update_by_x(x_plus_delta)
            ll_delta, _ = fn(edge_derivative = False)
            d_estimate = (ll_delta - ll) / delta           
            other_derivs.append(d_estimate)
            # restore self.x
            self.update_by_x(x)
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

    def objective_and_gradient(self, display, x):
        self.update_by_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        return f, g

    def Clock_wrap(self, display, x_clock):
        assert(self.clock)
        self.update_by_x_clock(x_clock)

        f, g = self.loglikelihood_and_gradient()
        
        # Now need to calculate the derivatives
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        other_derives, edge_derives = g[:-nEdge], g[-nEdge:]
        edge_to_derives = {self.edge_list[i] : edge_derives[i] for i in range(len(self.edge_list))}

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order

        Lr_derives = []  # used to store derivatives for the clock parameters L, r0, r1, ...
        Lr_derives.append(sum(edge_derives))  # dLL/dL = sum(all derives)
        Lr_derives.append(edge_to_derives[out_group_branch] * 2 / (self.x_Lr[1] - 2)
                          + sum(edge_derives))

        for i in range(2, len(self.x_Lr)):  # r(i-1)
            Lr_derives.append( edge_to_derives[('N' + str(i - 2), 'N' + str(i - 1))] * self.x_Lr[i] / (self.x_Lr[i] - 1) # 
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


    def get_mle(self, display = True):
        bnds = [(None, -0.05)] * 3
        if not self.clock:
            self.update_by_x()
            f = partial(self.objective_and_gradient, display)
            guess_x = self.x            
            bnds.extend([(None, None)] * (len(self.x) - 3))
            
        else:
            self.update_by_x_clock()  # TODO: change force for blen in x_clock
            f = partial(self.Clock_wrap, display)
            guess_x = self.x_clock
            bnds.extend([(None, None)] * (len(self.x_clock) - 2 - (len(self.edge_to_blen) / 2 + 1)))
            bnds.extend([(-10, 0.0)] * (len(self.edge_to_blen) / 2))        

        result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        print (result)
        return result

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
            ExpectedGeneconv = {self.edge_list[i] : j_out['responses'][0][i] for i in range(len(self.edge_list))}
            return ExpectedGeneconv
        else:
            print 'Need to implement this for old package'

    def get_ExpectedNumGeneconv(self):
        self.ExpectedGeneconv = self._ExpectedNumGeneconv()
    
    def save_to_file(self, file_name = None, path = './'):
        if file_name == None:
            file_name = self.Model + '_' + '_'.join(self.paralog)
            if self.clock:
                file_name += '_clock.p'
            else:
                file_name += '_nonclock.p'
        save_file = path + file_name
        self.ll = self._loglikelihood()[0]
        print 'x = ', self.x, 'x_clock = ', self.x_clock
        save_info = dict(
            Model = self.Model,
            x     = self.x,
            x_clock = self.x_clock,
            edge_to_blen = self.edge_to_blen,
            edge_list = self.edge_list,
            ExpectedGeneconv = self.ExpectedGeneconv,
            ll    = self.ll,
            newicktree = self.newicktree,
            alignment_file = self.seqloc,
            paralog = self.paralog,
            clock = self.clock,
            Force = self.Force
            )
        pickle.dump(save_info, open(save_file, 'wb+'))  # use pickle to save the class which can be easily reconstructed by pickle.load()

def main(args):
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    path = './NewPackageNewRun/'

    print 'Now calculate MLE for pair', paralog
    x = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574,
        1.19617746, -3.90910406, -3.87159909, -4.30005794, -6.28014701,
       -2.83555499, -3.29056063, -2.92359461, -3.49891453, -3.97438916,
       -3.10701318, -4.59465356, -2.01786889])
    
    if hasattr(args, 'Force'):
        Force = args.Force
        if Force != None:
            path += 'Force_'
    else:
        Force = None

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False)
    test.update_by_x(x)
    
    result = test.get_mle(display = True)
    test.get_ExpectedNumGeneconv()
    test.save_to_file(path = path)

    x_clock = np.array([-0.64430238, -0.48627304, -0.93852343, 1.76793238, -2.45061608, 1.23836258,
                  -2.11088631, -0.15313429, -0.23202426, -0.63091448, -0.25441366, -0.20629961,-0.23301088])
    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = True)
    test2.update_by_x_clock(x_clock)
    result = test2.get_mle(display = True)
    test2.get_ExpectedNumGeneconv()
    test2.save_to_file(path = path)
    

if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--Force', type = ast.literal_eval, help = 'Parameter constraints')
##    
##    main(parser.parse_args())

    paralog1 = 'YML026C'
    paralog2 = 'YDR450W'
    #paralog1 = 'ECP'
    #paralog2 = 'EDN'

    Force    = {5:0.0}

    paralog = [paralog1, paralog2]
    alignment_file = '../PairsAlignemt/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    #alignment_file = '../data/cleaned_input_data.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    #newicktree = '../data/input_tree.newick'
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False)#, nnsites = 5)
    #test.get_mle(display = False)
    #test.get_ExpectedNumGeneconv()
    #print test.ExpectedGeneconv
    x = np.array([-0.64353555, -0.48264002, -0.93307917,  1.76977596, -2.44028574, 1.19617746,
                   -3.90910406, -2.01786889,  # N0
                   -2.92359461, -2.83555499,  # N1
                   -4.30005794, -3.29056063,  # N2
                   -4.59465356, -3.49891453,  # N3
                   -6.28014701, -3.10701318,  # N4
                   -3.87159909, -3.97438916]) # N5
                 
    test.update_by_x(x = x)

    print test._loglikelihood()  # should be -1513.0033676643861


##    paralog1 = 'YDR502C'
##    paralog2 = 'YLR180W'
##
##
##    Force    = {5:0.0}
##
##    paralog = [paralog1, paralog2]
##    alignment_file = '../PairsAlignemt/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##
##    newicktree = '../PairsAlignemt/YeastTree.newick'
##    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False)#, nnsites = 5)
##
##
##    x = np.array([-0.76374179, -0.56512517, -0.8661595 ,  1.33286957, -3.12679708, -0.55435137,
##                  -1.59080985, -1.35214307,   # N0
##                  -1.48481827, -0.5117971 ,   # N1
##                  -2.4550277 , -2.02888192,   # N2
##                  -2.47147681, -1.7627708 ,   # N3
##                  -2.62667816, -1.64821228,   # N4
##                  -2.38731535, -2.41097359])  # N5
##                  
##    test2.update_by_x(x = x)
##
##    print test2._loglikelihood()
##        
        
        
