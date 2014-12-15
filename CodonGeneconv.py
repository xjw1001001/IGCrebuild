# Uses Alex Griffing's JsonCTMCTree package for likelihood and gradient calculation
# Re-written of my previous CondonBased2Repeats class
from CodonGeneconFunc import *

class CodonGeneconv:
    def __init__(self, tree_newick, alignment, paralog, Model = 'MG94', nnsites = None):
        self.newicktree = tree_newick  # newick tree file loc
        self.seqloc = alignment # multiple sequence alignment, now need to remove gap before-hand
        self.paralog = paralog  # paralog list
        self.nsites = nnsites  # number of sites in the alignment used for calculation
        self.Model = Model


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
         # self.pi         real values
         # self.kappa      real values
         # self.omega      real values
         # self.tau        real values
        self.processes = self.get_MG94Geneconv_and_MG94()
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
        x_process = np.log(x_process)        
        setattr(CodonGeneconv, 'x_process', x_process)
        
        x_rates = np.array([ self.edge_to_blen[edge] for edge in self.edge_to_blen.keys()])
        x_rates = np.log(x_rates)
        setattr(CodonGeneconv, 'x_rates', x_rates)

        x = np.concatenate((x_process, x_rates))
        setattr(CodonGeneconv, 'x', x)

    def update_by_x(self, x = None):
        k = len(self.edge_to_blen)
        if x != None:
            self.x = x
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        self.unpack_x_process()
        self.unpack_x_rates()
        

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
        x_process = np.exp(self.x_process)
        pi_a = x_process[0] * x_process[1]
        pi_c = (1 - x_process[0]) * x_process[2]
        pi_g = x_process[0] * (1 - x_process[1])
        pi_t = (1 - x_process[0]) * (1 - x_process[2])
        self.pi = [pi_a, pi_c, pi_g, pi_t]
        self.kappa = x_process[3]
        self.omega = x_process[4]
        self.tau = x_process[5]

    def unpack_x_rates(self):
        x_rates = np.exp(self.x_rates)
        assert(len(x_rates) == len(self.edge_to_blen))
        for it in range(len(self.edge_to_blen)):
            edge = self.edge_to_blen.keys()[it]
            self.edge_to_blen[edge] = x_rates[it]
        

    def nts_to_codons(self):
        for name in self.name_to_seq.keys():
            assert(len(self.name_to_seq[name]) % 3 == 0)
            tmp_seq = [self.name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(self.name_to_seq[name]) / 3 )]
            self.name_to_seq[name] = tmp_seq
            
        
    def get_data(self, Model):
        seq_dict = SeqIO.to_dict(SeqIO.parse( self.seqloc, "fasta" ))
        setattr(CodonGeneconv, 'name_to_seq', {name:str(seq_dict[name].seq) for name in seq_dict.keys()})

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
                observation = obs_to_state[self.name_to_seq[name][site]]
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

    def vec_get_MG94Geneconv(self):
        pair_list = [a + b for a, b in product(self.codon_nonstop, repeat = 2)]
        pair_from_matrix = np.repeat(pair_list, len(pair_list)).reshape((len(pair_list), len(pair_list)))
        process_geneconv = vec_get_MG94GeneconvRate(pair_from = pair_from_matrix, pair_to = pair_list, pi = self.pi,
                                   kappa = self.kappa, omega = self.omega, codon_table = self.codon_table,
                                   tau = self.tau, codon_to_state = self.codon_to_state)
    
    def get_HKYGeneconv(self, x = None):
        row = []
        col = []
        rate_geneconv = []
        rate_basic = []

        self.update_by_x(x)

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

        self.update_by_x(x)

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
        k = len(self.edge_to_blen)
        self.x_process, self.x_rates = self.x[:-k], self.x[-k:]
        if self.Model == 'MG94':
            self.processes = self.get_MG94Geneconv_and_MG94()
        elif self.Model == 'HKY':
            self.processes = self.get_HKYGeneconv()

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

        status = j_ll['status']
        feasibility = j_ll['feasibility']

        if status != 'success' or not feasibility:
            print 'results:'
            print j_ll
            print
            raise Exception('Encountered some problem in the calculation of log likelihood and its derivatives')

        ll, edge_derivs = j_ll['log_likelihood'], j_ll['edge_derivatives']

        return ll, edge_derivs
        
        
    def loglikelihood_and_gradient(self, display = False):
        '''
        Modified from Alex's objective_and_gradient function in ctmcaas/adv-log-likelihoods/mle_geneconv_common.py
        '''
        delta = 1e-8
        x = deepcopy(self.x)  # store the current x array

        ll, edge_derivs = self._loglikelihood(edge_derivative = True)
        
        m = len(self.x) - len(self.edge_to_blen)

        # use finite differences to estimate derivatives with respect to these parameters
        other_derivs = []
        
        for i in range(m):
        #for i in range(m-1, -1, -1):
            x_plus_delta = np.array(x)
            x_plus_delta[i] += delta
            self.update_by_x(x_plus_delta)
            ll_delta, _ = self._loglikelihood(edge_derivative = False)
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

        f = -ll
        g = -np.concatenate((other_derivs, edge_derivs))
        return f, g
        #ll, edge_derivs =
    def objective_and_gradient(self, display, x):
        self.update_by_x(x)
        f, g = self.loglikelihood_and_gradient(display = display)
        return f, g

    def get_mle(self, display = True):
        guess_x = self.x
        bnds = [(None, -0.05)] * 3
        bnds.extend([(None, None)] * (len(self.x) - 3))
        f = partial(self.objective_and_gradient, display)
        result = scipy.optimize.minimize(f, guess_x, jac = True, method = 'L-BFGS-B', bounds = bnds)
        print (result)
        return result

        


if __name__ == '__main__':
    #alignment_file = '/Users/xji3/Genconv/PairsAlignemt/YDR502C_YLR180W/YDR502C_YLR180W_input.fasta'
    alignment_file = '/Users/xji3/Genconv/PairsAlignemt/YAL056W_YOR371C/YAL056W_YOR371C_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    paralog = ['YAL056W', 'YOR371C']
    test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', nnsites = 10)
    #test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', nnsites = 10)
##    x = np.array([-0.4671754 , -0.47657272, -1.05338269,  0.58005615, -1.62906811,
##       -1.42286024,  0.3567071 , -0.85263062, -0.38357827, -0.65378769,
##       -0.51712782, -0.33653191,  0.60581195,  0.84021518, -0.65989123,
##       -0.62867968, -0.55471522,  0.3931776])
##
##    x = np.array([-0.52856124, -0.50895601, -1.08978484,  0.14326525, -2.45656832,
##       -1.62703482,  0.47515258, -1.24235701, -0.62424348, -1.00954971,
##       -0.80980144, -0.54295938,  0.86255682,  1.21904872, -0.95233185,
##       -0.96061732, -0.88000097,  0.57393683])
##
##    x = np.array([-0.33192463, -0.38047998, -0.945599  ,  0.09092072, -1.54660303,
##       -1.56649691,  0.70287302, -1.50665575, -0.86620656, -1.41793859,
##       -0.05129703, -0.71300579,  1.82643417,  1.03687331, -1.40222114,
##       -1.32922251, -1.22190485,  1.07804555])
##    test.update_by_x(x)


    print 'Now calculate likelihood'
    #j_ll = test.loglikelihood_and_gradient()

    result = test.get_mle(display = True)
    #result = test.get_mle(display = True)

    #cProfile.run('result = test.get_mle()')

##    print get_MG94BasicRate('ATG','AAG',test.pi,test.kappa, test.omega, test.codon_table)
##    print vec_get_MG94BasicRate(ca = ['ATG','TAG'],cb = ['AAG','ATG'],pi = test.pi,kappa = test.kappa, omega = test.omega, codon_table = test.codon_table)
##
##    print get_MG94GeneconvRate('ATGAAG','ATGATG', test.pi, test.kappa, test.omega, test.codon_table, test.tau, test.codon_to_state)
##    print vec_get_MG94GeneconvRate(pair_from = ['ATGAAG'], pair_to = ['ATGATG', 'ATGAAT'], pi = test.pi,
##                                   kappa = test.kappa, omega = test.omega, codon_table = test.codon_table,
##                                   tau = test.tau, codon_to_state = test.codon_to_state)
