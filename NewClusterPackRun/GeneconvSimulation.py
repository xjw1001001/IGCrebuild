from CodonGeneconFunc import *
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


        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.state_to_codon = {i:a.upper() for (i, a) in enumerate(self.codon_nonstop)}
        self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}
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

        # Prior distribution on the root
        self.prior_feasible_states  = None
        self.prior_distribution     = None

        # Node sequence
        self.node_to_sequence = None

        self.initiate()

    def initiate(self):
        self.get_tree()
        self.unpack_x()

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
                    Tgeneconv12 = self.tau[0] * self.omega
                    Tgeneconv21 = self.tau[1] * self.omega
                else:
                    Tgeneconv12 = self.tau[0]
                    Tgeneconv21 = self.tau[1]
                rate_geneconv.append(Qb + Tgeneconv12 * self.get_GC_fitness(cb, ca))
                rate_basic.append(0.0)
                
                # (ca, cb) to (cb, cb)
                row.append((sa, sb))
                col.append((sb, sb))
                Qb = Qbasic[sa, sb]
                rate_geneconv.append(Qb + Tgeneconv21 * self.get_GC_fitness(ca, cb))
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

        mat_row = [self.pair_to_state[(self.state_to_codon[i[0]], self.state_to_codon[i[1]])] for i in row]
        mat_col = [self.pair_to_state[(self.state_to_codon[i[0]], self.state_to_codon[i[1]])] for i in col]
        Geneconv_mat = csr_matrix((process_geneconv['rate'], (mat_row, mat_col)), shape = (61**2, 61**2), dtype = float)
        basic_mat = csr_matrix((process_basic['rate'], (mat_row, mat_col)), shape = (61**2, 61**2), dtype = float)
        return Geneconv_mat, basic_mat

    def draw_from_distribution(self, prob, size):
        if self.Model == 'MG94':
            values = self.pair_to_state.keys()
        bins = np.add.accumulate(prob)
        return values[np.digitize(np.random.random_sample(size), bins)]


if __name__ == '__main__':
    
    paralog1 = 'YLR406C'
    paralog2 = 'YDL075W'

    paralog = [paralog1, paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    log_omega = np.log(8.562946237878546474e-02)
    x = np.exp([-0.69356832, -0.50733147, -0.88584866,  2.33732053, log_omega,  2.08070563,
       -3.17876664, -3.06491686, -3.3344096 , -3.59554012, -5.05790222,
       -3.70204719, -5.24464693, -4.11413797, -6.24278987, -3.75290647,
       -4.5716839 , -4.44936138])
    
    test = SimGeneconv( newicktree, paralog, x, Model = 'MG94', nnsites = 2)

    self = test

    g, b = self.get_MG94Geneconv_and_MG94()
    prob = np.array(g.getrow(0).todense())[0,:]
    self.draw_from_distribution(prob, 1)
