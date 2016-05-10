from IGCexpansion.CodonGeneconFunc import *
from scipy.sparse import csr_matrix


class SimGeneconv:
    def __init__(self, tree_newick, paralog, x, Model = 'MG94', nnsites = 1000, Dir = False, gBGC = False):
        self.newicktree = tree_newick
        self.paralog      = paralog
        self.Model        = Model
        self.nsites       = nnsites
        self.Dir          = Dir
        self.gBGC         = gBGC
        
        
        self.tree         = None
        self.edge_to_blen = None
        self.node_to_num  = None
        self.num_to_node  = None
        self.edge_list    = None

        self.outgroup     = [('N0', 'kluyveri')]
        self.leaves       = None


        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.state_to_nt    = {i:a for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.state_to_codon = {i:a.upper() for (i, a) in enumerate(self.codon_nonstop)}
        if self.Model == 'MG94':
            self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}
        elif self.Model == 'HKY':
            self.pair_to_state  = {pair:i for i, pair in enumerate(product('ACGT', repeat = 2))}
        self.state_to_pair  = {self.pair_to_state[pair]:pair for pair in self.pair_to_state}


        # Rate matrix related variable
        self.x_process      = None      # values of process parameters 
        self.x_rates        = None      # values of blen 
        self.x              = x         # x_process + x_rates
        self.pi             = None      # real values
        self.kappa          = None      # real values
        self.omega          = None      # real values
        self.tau            = None      # real values  Tau12, Tau21
        self.gamma          = 0.0
        self.L              = 200.0       # average gene conversion length

        # Prior distribution on the root
        self.prior_feasible_states  = None
        self.prior_distribution     = None

        # Node sequence
        self.node_to_sequence = None
        self.node_to_sim      = {}

        # Rate Matrices
        self.Geneconv_mat     = None
        self.Basic_mat        = None
        self.IGC_mat          = None

        self.initiate()

    def initiate(self):
        self.get_tree()
        self.unpack_x()
        if self.Model == 'MG94':
            self.Geneconv_mat, self.Basic_mat, self.IGC_mat = self.get_MG94Geneconv_and_MG94()
        elif self.Model == 'HKY':
            self.Geneconv_mat, self.Basic_mat = self.get_HKYGeneconv_and_HKY()

    def unpack_x(self):
        k = len(self.edge_to_blen)
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        self.unpack_x_process()
        self.unpack_x_rates()

    def unpack_x_rates(self):  # TODO: Change it to fit general tree structure rather than cherry tree
        x_rates = self.x_rates
        
        assert(len(x_rates) == len(self.edge_to_blen))

        for edge_it in range(len(self.edge_list)):
            self.edge_to_blen[self.edge_list[edge_it]] = x_rates[edge_it]

    def unpack_x_process(self):
        x_process = self.x_process
        # %AG, %A, %C, kappa, tau
        num_of_para = 5
        if self.Dir:
            num_of_para += 1
        if self.gBGC:
            num_of_para += 1
        if self.Model == 'MG94':
            num_of_para += 1
        assert(len(x_process) == num_of_para)

        if self.Model == 'MG94':
            # x_process[] = %AG, %A, %C, kappa, omega, tau
            
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            self.omega = x_process[4]
            if self.Dir:
                self.tau = x_process[5:7]
            else:
                self.tau = [x_process[5]] * 2
            if self.gBGC:
                self.gamma = self.x_process[-1]

        elif self.Model == 'HKY':
            # x_process[] = %AG, %A, %C, kappa, tau
            pi_a = x_process[0] * x_process[1]
            pi_c = (1 - x_process[0]) * x_process[2]
            pi_g = x_process[0] * (1 - x_process[1])
            pi_t = (1 - x_process[0]) * (1 - x_process[2])
            self.pi = [pi_a, pi_c, pi_g, pi_t]
            self.kappa = x_process[3]
            if self.Dir:
                self.tau = x_process[4:6]
            else:
                self.tau = [x_process[4]] * 2
            if self.gBGC:
                self.gamma = self.x_process[-1]            

        self.node_to_sequence = {node:[] for node in self.node_to_num.keys()}
        self.get_prior()
            
            
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
        self.leaves = list(leaves)
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
        #self.sim_root()

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

    def get_GC_fitness(self, ca, cb):
        assert(len(ca) == len(cb)) # length should be the same
        # ca got copied by cb. 
        k = cb.count('G') + cb.count('C') - ca.count('G') - ca.count('C')
        if self.gamma == 0.0:
            return 1.0
        else:
            if k == 0:
                return 1.0
            else:
                return (k * self.gamma) / (1 - np.exp(- k * self.gamma))
        
    def get_MG94Geneconv_and_MG94(self):
        Qbasic = self.get_MG94Basic()
        row = []
        col = []
        rate_geneconv = []
        rate_basic    = []
        rate_IGC      = []

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
                        rate_IGC.append(0.0)

                    # (ca, cb) to (cc, cb)
                    Qb = Qbasic[sa, sc]
                    if Qb != 0:
                        row.append((sa, sb))
                        col.append((sc, sb))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)
                        rate_IGC.append(0.0)

                        
                # (ca, cb) to (ca, ca)
                row.append((sa, sb))
                col.append((sa, sa))
                Qb = Qbasic[sb, sa]
                if isNonsynonymous(cb, ca, self.codon_table):
                    Tgeneconv12 = self.tau[0] * self.omega
                    Tgeneconv21 = self.tau[1] * self.omega
                else:
                    Tgeneconv12 = self.tau[0]
                    Tgeneconv21 = self.tau[1]
                rate_geneconv.append(Qb + Tgeneconv12 * self.get_GC_fitness(cb, ca))
                rate_basic.append(0.0)
                rate_IGC.append(Tgeneconv12 * self.get_GC_fitness(cb, ca))
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qb = Qbasic[sa, sb]
                rate_geneconv.append(Qb + Tgeneconv21 * self.get_GC_fitness(ca, cb))
                rate_basic.append(0.0)
                rate_IGC.append(Tgeneconv21 * self.get_GC_fitness(ca, cb))

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
                        rate_IGC.append(0.0)
                        
                    # (ca, ca) to (cc, ca)
                        row.append((sa, sb))
                        col.append((sc, sa))
                        rate_geneconv.append(Qb)
                        rate_basic.append(0.0)
                        rate_IGC.append(0.0)

                    # (ca, ca) to (cc, cc)
                        row.append((sa, sb))
                        col.append((sc, sc))
                        rate_geneconv.append(0.0)
                        rate_basic.append(Qb)
                        rate_IGC.append(0.0)
                
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
        process_IGC = dict(
            row = row,
            col = col,
            rate = rate_IGC
            )

        mat_row = [self.pair_to_state[(self.state_to_codon[i[0]], self.state_to_codon[i[1]])] for i in row]
        mat_col = [self.pair_to_state[(self.state_to_codon[i[0]], self.state_to_codon[i[1]])] for i in col]
        Geneconv_mat = csr_matrix((process_geneconv['rate'], (mat_row, mat_col)), shape = (61**2, 61**2), dtype = float)
        basic_mat = csr_matrix((process_basic['rate'], (mat_row, mat_col)), shape = (61**2, 61**2), dtype = float)
        IGC_mat = csr_matrix((process_IGC['rate'], (mat_row, mat_col)), shape = (61**2, 61**2), dtype = float)
        return Geneconv_mat, basic_mat, IGC_mat

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
        

    def get_HKYGeneconv_and_HKY(self):
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

                if (na != nc and nb!= nd) or (na == nc and nb == nd):
                    GeneconvRate = 0.0
                if na ==nc and nb != nd:
                    Qb = Qbasic['ACGT'.index(nb), 'ACGT'.index(nd)]
                    if na == nd:
                        GeneconvRate = Qb + self.tau[0] * self.get_GC_fitness(nb, na)
                    else:
                        GeneconvRate = Qb
                if nb == nd and na != nc:
                    Qb = Qbasic['ACGT'.index(na), 'ACGT'.index(nc)]
                    if nb == nc:
                        GeneconvRate = Qb + self.tau[1] * self.get_GC_fitness(na, nb)
                    else:
                        GeneconvRate = Qb

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
        mat_row = [self.pair_to_state[(self.state_to_nt[i[0]], self.state_to_nt[i[1]])] for i in row]
        mat_col = [self.pair_to_state[(self.state_to_nt[i[0]], self.state_to_nt[i[1]])] for i in col]
        Geneconv_mat = csr_matrix((process_geneconv['rate'], (mat_row, mat_col)), shape = (4**2, 4**2), dtype = float)
        basic_mat = csr_matrix((process_basic['rate'], (mat_row, mat_col)), shape = (4**2, 4**2), dtype = float)
        return Geneconv_mat, basic_mat


    def draw_from_distribution(self, prob, size,values = None, prob_IGC = None):
        if values == None:
            values = [self.state_to_pair[i] for i in range(len(self.state_to_pair))]
        bins = np.add.accumulate(prob)
        out_state = None
        num_IGC = 0
        if size == 1:
            event = np.digitize(np.random.random_sample(size), bins)
            if prob_IGC is not None:
                #print '', event, 'prob dist = ', prob_IGC[event], prob[event], self.tau[0]
                if prob_IGC[event] > 0.0:
                    if (prob_IGC[event] / prob[event]) > 1.0:
                        assert( 0 ==1)
                    if np.random.uniform() < (prob_IGC[event] / prob[event]):
                        num_IGC += 1
                        
            if prob_IGC != None:
                return values[event], num_IGC
            else:
                return values[event]
        else:
            event_list = np.digitize(np.random.random_sample(size), bins)
            for event in event_list:
                if prob_IGC != None:
                    print prob_IGC
                    if prob_IGC[event] > 0.0:
                        if np.random.uniform() < (prob_IGC[event] / prob[event]):
                            num_IGC += 1
            if prob_IGC != None:
                return [values[i] for i in event_list], num_IGC
            else:
                return [values[i] for i in event_list]
                

    def sim_one_branch(self, starting_state, rate_matrix, IGC_matrix, blen):
        num_IGC = 0
        num_All = 0
        cummulate_time = 0.0
        #starting_state = self.pair_to_state[starting_pair]
        prob = np.array(rate_matrix.getrow(starting_state).todense())[0,:]
        prob_IGC = np.array(IGC_matrix.getrow(starting_state).todense())[0,:]
        while(cummulate_time < blen):
            #print 'cummulated time = ', cummulate_time, 'starting State = ', starting_state, 'blen = ', blen
            if sum(prob) == 0.0:
                break
            cummulate_time += np.random.exponential(1.0 / sum(prob))
            #print 'Next event time = ', cummulate_time
            if cummulate_time > blen:
                break
            else:
                starting_state, add_IGC_num = self.draw_from_distribution(prob = deepcopy(prob) / sum(prob), size = 1, values = range(len(prob)), prob_IGC = deepcopy(prob_IGC)/sum(prob))
                num_IGC += add_IGC_num
                num_All += 1
                prob = np.array(rate_matrix.getrow(starting_state).todense())[0,:]
                prob_IGC = np.array(IGC_matrix.getrow(starting_state).todense())[0,:]
                #print prob
                #print prob_IGC
                #print sum(prob_IGC), sum(prob), num_All, num_IGC
                #starting_state = self.pair_to_state[starting_pair]

        return starting_state, num_IGC, num_All #starting_pair#, starting_state, cummulate_time


    def sim_root(self):
        if self.Model == 'MG94':
            seq = self.draw_from_distribution(self.prior_distribution, self.nsites, self.codon_nonstop)
        elif self.Model == 'HKY':
            seq = self.draw_from_distribution(self.prior_distribution, self.nsites, 'ACGT')
        self.node_to_sequence['N0'] = np.array([self.pair_to_state[(i, i)] for i in seq])
        self.node_to_sim['N0'] = [self.node_to_sequence['N0'], 0, 0]

    #vec_sim_one_branch = np.vectorize(self.sim_one_branch, excluded = ['rate_matrix', 'blen'])

    def sim(self):
        self.sim_root()
        for edge in self.edge_list:
            print edge
            if edge in self.outgroup:
                rate_mat = self.Basic_mat
            else:
                rate_mat = self.Geneconv_mat

            IGC_mat = self.IGC_mat

            end_seq = []
            num_IGC = 0
            num_All = 0
            for site in self.node_to_sequence[edge[0]]:
                #print site
                site_seq, add_num_IGC, add_num_All = self.sim_one_branch(site, rate_mat, IGC_mat, self.edge_to_blen[edge])
                num_IGC += add_num_IGC
                num_All += add_num_All
                end_seq.append(site_seq)
            self.node_to_sim[edge[1]] = [end_seq, num_IGC, num_All]
            self.node_to_sequence[edge[1]] = end_seq

        for node in self.node_to_sim.keys():
            seq1 = ''.join([self.state_to_pair[i][0] for i in self.node_to_sim[node][0]])
            seq2 = ''.join([self.state_to_pair[i][1] for i in self.node_to_sim[node][0]])
            self.node_to_sequence[node] = (seq1, seq2)

    def sim_tract(self):
        self.sim_root()
        self.node_to_sequence['N0'] = [''.join([self.state_to_pair[i][0] for i in self.node_to_sequence['N0']]),
                                       ''.join([self.state_to_pair[i][1] for i in self.node_to_sequence['N0']])]
        for edge in self.edge_list:
            print edge
            if edge in self.outgroup:
                rate_mat = self.Basic_mat
            else:
                rate_mat = self.Geneconv_mat

            end_seq = []
            for site in self.node_to_sequence[edge[0]]:
                #print site
                end_seq.append(self.sim_one_branch(site, rate_mat, self.edge_to_blen[edge]))
            self.node_to_sequence[edge[1]] = end_seq

        for node in self.node_to_sequence.keys():
            seq1 = ''.join([self.state_to_pair[i][0] for i in self.node_to_sequence[node]])
            seq2 = ''.join([self.state_to_pair[i][1] for i in self.node_to_sequence[node]])
            self.node_to_sequence[node] = (seq1, seq2)        

    def sim_one_branch_tract(self, starting_seq_pair, blen):
        assert(len(starting_seq_pair[0]) == len(starting_seq_pair[1]))
        # TODO: output total number of geneconv events
        cummulate_time = 0.0
        Qbasic = self.get_MG94Basic()
        diagonal_rates = Qbasic.sum(axis = 1)

        # prepare rate array for sampling exponential waiting time
        seq_rate_1 = [diagonal_rates[self.codon_to_state[starting_seq_pair[0][i * 3:(i + 1) * 3]]] for i in range(len(starting_seq_pair[0]) / 3)]
        seq_rate_2 = [diagonal_rates[self.codon_to_state[starting_seq_pair[1][i * 3:(i + 1) * 3]]] for i in range(len(starting_seq_pair[1]) / 3)]

        # Now prepare "rate matrix" for the gene conversion event
        # row number is the starting position, column number is the ending position
        S = len(starting_seq_pair[0])
        J = np.zeros((S, S))  # copying my formula part from Alex's script
        l = np.arange(S+1, dtype=float)
        p = (1/self.L) * (1 - 1/self.L)**l  # geometric distribution
        for i in range(S):
            for j in range(i, S):
                l = j - i  # This is the number of failures in Geometric distribution
                # start and end inside the gene with length l <= S - 2
                if 0 < i <= j < S - 1:
                    J[i, j] += p[l]
                # start before gene, end within gene with length l <= S - 1
                if 0 == i <= j < S - 1:
                    J[i, j] += p[l] * self.L
                # start within gene, end after gene with length l <= S - 1
                if 0 < i <= j == S - 1:
                    J[i, j] += p[l] * self.L
                # start before gene, end after gene with length l = S
                if 0 == i <= j == S - 1:
                    J[i, j] += p[l] * self.L ** 2                
        # Now assign initiation rate using tau parameter
        init_rate = sum(self.tau) / 2.0 / self.L
        J *= init_rate

        mut_prob = np.array(seq_rate_1 + seq_rate_2)
        IGC_prob = np.ravel(J)
        total_rate = mut_prob.sum() + IGC_prob.sum() * 2
        
        while(cummulate_time < blen):
            #Now sample exponential distributed waiting time for next event (either point mutation or IGC)    
            if total_rate == 0.0:
                break
            cummulate_time += np.random.exponential(1.0 / total_rate)
            if cummulate_time > blen:
                break
            else:
                # an event happening
                event = self.draw_from_distribution(np.concatenate((mut_prob, IGC_prob, IGC_prob)) / total_rate, 1, range(len(mut_prob) + len(IGC_prob) * 2))
                #event = np.floor(np.random.uniform() * (len(mut_prob) + len(IGC_prob) * 2))
                if event < len(mut_prob):  # point mutation event
                    if event < len(seq_rate_1):  # mutation on seq 1
                        prob = np.array(Qbasic.getrow(self.codon_to_state[starting_seq_pair[0][event * 3:(event + 1) * 3]])[0,:])
                        new_codon = self.codon_nonstop[self.draw_from_distribution(deepcopy(prob) / sum(prob), 1, range(len(prob)))]
                        starting_seq_pair[0] = starting_seq_pair[0][:event * 3] + new_codon + starting_seq_pair[0][(event + 1) * 3:]
                    else:
                        event = event - len(seq_rate_1)
                        prob = np.array(Qbasic[self.codon_to_state[starting_seq_pair[1][event * 3:(event + 1) * 3]],:])
                        new_codon = self.codon_nonstop[self.draw_from_distribution(deepcopy(prob) / sum(prob), 1, range(len(prob)))]
                        starting_seq_pair[1] = starting_seq_pair[1][:event * 3] + new_codon + starting_seq_pair[1][(event + 1) * 3:]
                        
                else:  # IGC event
                    event -= len(mut_prob)
                    if event > len(IGC_prob):  # seq 2 being donor seq
                        event -= len(IGC_prob)
                        init_pos = int(np.floor((event / S)))  # starting position
                        stop_pos = event - init_pos * S  # stop position
                        # now copy from donor seq to recipient seq
                        starting_seq_pair[0] = starting_seq_pair[0][:init_pos] + starting_seq_pair[0][init_pos:(stop_pos + 1)] + starting_seq_pair[0][(stop_pos + 1):]
                    else:
                        init_pos = int(np.floor((event / S)))  # starting position
                        stop_pos = event - init_pos * S  # stop position
                        # now copy from donor seq to recipient seq
                        starting_seq_pair[1] = starting_seq_pair[1][:init_pos] + starting_seq_pair[1][init_pos:(stop_pos + 1)] + starting_seq_pair[1][(stop_pos + 1):]                        
                    
        return starting_seq_pair#starting_pair#, starting_state, cummulate_time
                
    def output_seq(self, path = './simulation/', sim_num = None, Force = False, clock = False):
        file_name = '_'.join(self.paralog) + '_' + self.Model
        if self.Dir:
            file_name += '_dir'
        if self.gBGC:
            file_name += '_gBGC'
        if clock:
            file_name += '_clock'
        else:
            file_name += '_nonclock'
        if Force:
            file_name += '_force'

        if sim_num != None:
            file_name += '_Sim_' + str(sim_num)

        output_file = path + file_name + '.fasta'
        outgroup_leaves = [i[1] for i in self.outgroup]
        with open(output_file, 'w+') as f:
            for node in self.leaves:
                if node in outgroup_leaves:
                    f.write('>' + node + self.paralog[0] + '\n')
                    f.write(self.node_to_sequence[node][0] + '\n')
                else:
                    f.write('>' + node + self.paralog[0] + '\n')
                    f.write(self.node_to_sequence[node][0] + '\n')
                    f.write('>' + node + self.paralog[1] + '\n')
                    f.write(self.node_to_sequence[node][1] + '\n')
                    
                
   


if __name__ == '__main__':
    
    paralog1 = 'YLR406C'
    paralog2 = 'YDL075W'

##    paralog1 = 'YML026C'
##    paralog2 = 'YDR450W'

    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'

    x = np.exp(np.loadtxt(open('./save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt', 'r')))

    test = SimGeneconv(newicktree, paralog, x, Model = 'MG94', nnsites = 400, Dir = False, gBGC = False)
    self = test
##    self.sim_root()
##    self.tau = [0.0] * 2
##    IGC_mat = self.Geneconv_mat - self.get_MG94Geneconv_and_MG94()[0]
##    edge = self.edge_list[0]
##    site = self.node_to_sequence[edge[0]][0]
##    rate_matrix = self.Geneconv_mat
##    a,b,c = test.sim_one_branch(site, rate_matrix, IGC_mat, 1.0)
    test.sim()
    #test.sim_root()
    test.sim_one_branch(test.node_to_sim['N0'][0][0], test.Geneconv_mat, test.IGC_mat, 1.0)
    print test.pi
    print (test.node_to_sequence['N0'][0].count('A') + 0.0) / (test.nsites * 3.0), (test.node_to_sequence['N0'][0].count('C') + 0.0) / (test.nsites * 3.0), (test.node_to_sequence['N0'][0].count('G') + 0.0) / (test.nsites * 3.0), (test.node_to_sequence['N0'][0].count('T') + 0.0) / (test.nsites * 3.0)
    print sum([test.node_to_sim[i][1] + 0.0 for i in test.node_to_sim if i != 'kluyveri']) / sum([test.node_to_sim[i][2] + 0.0 for i in test.node_to_sim if i != 'kluyveri'])

    test.output_seq(path = './Simulation/' + '_'.join(paralog) + '/', sim_num = sim_num)    
    #test = SimGeneconv( newicktree, paralog, x_gBGC_dir, Model = 'MG94', nnsites = 1000)
    #test = SimGeneconv( newicktree, paralog, x_HKY, Model = 'HKY', nnsites = 1000)

##    Dir = True
##    #Dir = False
##    gBGC = True
##    #gBGC = False
##    test = SimGeneconv( newicktree, paralog, x_gBGC_dir, Model = 'MG94', nnsites = 381, Dir = Dir, gBGC = gBGC)
##    self = test
##    self.sim_root()
##    self.node_to_sequence['N0'] = [''.join([self.state_to_pair[i][0] for i in self.node_to_sequence['N0']]),
##                                       ''.join([self.state_to_pair[i][1] for i in self.node_to_sequence['N0']])]
##    starting_seq_pair = self.node_to_sequence['N0']
##    blen = 0.1
    
##    test.sim()
##
##    print self.node_to_sequence['N0'][0].count('A'), self.node_to_sequence['N0'][0].count('A') / (self.nsites * 3.0)
##    print self.node_to_sequence['kluyveri'][0].count('A'), self.node_to_sequence['kluyveri'][0].count('A') / (self.nsites * 3.0)
##    print self.pi

    #self.output_seq()
