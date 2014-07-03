from CodonBased2Repeats import *
from numpy.random import random_sample, exponential, random

def weighted_values(values, probabilities, size=1):
    bins = np.add.accumulate(probabilities)
    return values[np.digitize(random_sample(size), bins)]

def test_simulation(Q,Tau_matrix,distn,state_list,state_to_num,t_end):
    s_list = []
    t_list = []
    s = weighted_values(state_list,distn)
    t = 0.0
    s_list.append(s)
    t_list.append(t)
    z = exponential(-Q[state_to_num[s],state_to_num[s]])
    geneconv_events = 0
    none_geneconv_events = 0
    while(t+z<t_end):
        t = t+z
        t_list.append(t)
        probabilities = Q[state_to_num[s],:]/(-Q[state_to_num[s],state_to_num[s]])
        probabilities[state_to_num[s]] =0.0
        s_to = weighted_values(state_list,probabilities)
        s_list.append(s_to)
        if not Tau_matrix[state_to_num[s],state_to_num[s_to]]==0:
            if random()<(Tau_matrix[state_to_num[s],state_to_num[s_to]]/Q[state_to_num[s],state_to_num[s_to]]):
                geneconv_events += 1
            else:
                none_geneconv_events += 1
        else:
            none_geneconv_events += 1
        s = s_to
        z = exponential(-Q[state_to_num[s],state_to_num[s]])
    return s_list,t_list, geneconv_events, none_geneconv_events
    
    
    


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

    Model = 'HKY'
    sim_Tao=1.0
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()

    guess = np.log([0.7,0.5,0.4,np.e**2])
    guess = np.append(guess,sim_Tao/2.0)

    guess_w_blen = [-1.0]*(len(test.edge_to_blen))
    guess_w_blen[0:6] = [0.0,-1.0,-2.0,-3.0,-4.0,-5.0]
    guess_w_blen.extend(guess)
    
    r1 = test.estimate(args, Model, guess_w_blen, est_blen = True)
##    test.modelnum = 3
##    test.Tao = 1.81968875933
##    test.para = np.array([ 0.4902551, 0.5770529, 0.50075712, 2.11381758])
##    test.edge_to_blen = {('N1', 'N2'): 0.0088094406225698342, ('N2', 'N3'): 0.010918921506121825, ('N3', 'Chimpanzee'): 0.0045642798967184643, ('N2', 'Orangutan'): 0.029862484609549191, ('N0', 'N1'): 0.070595334023490081, ('N3', 'Gorilla'): 0.0050102717128137648, ('N1', 'Macaque'): 0.051135575608256184, ('N0', 'Tamarin'): 0.10296832243092251}

    Tau_matrix = test.get_Tau_matrix()

    edge_to_M,edge_to_M_none = test.get_edge_to_M_Geneconv()

    expected_geneconv, expected_none_geneconv =  test.get_edge_to_expectednumchanges()
    print 'Expected number of Geneconv events per site'
    for v in expected_geneconv.keys():
        print v, expected_geneconv[v]

    print 'Expected number of none Geneconv events per site'
    for v in expected_none_geneconv.keys():
        print v, expected_none_geneconv[v]
        

    #a  = get_edge_to_distn2d(test.treetopo,test.get_P(test.edge_to_blen,test.modelnum,test.para,test.Tao), 'N0', d,

    self = test
    NumRepeats = 2
    Q_post_duplication,d = self.get_Q_and_distn(self.modelnum, self.para, self.Tao, NumRepeats)
    expected_rate = self.get_expected_rate(Q_post_duplication)
    normalizing_factor = NumRepeats/expected_rate
    Q_post_duplication_normalized = normalizing_factor*Q_post_duplication
    
    Q_pre_duplication,d_pre = self.get_Q_and_distn(self.modelnum, self.para, self.Tao, 1)
    Tau_matrix = self.get_Tau_matrix(normalizing_factor,NumRepeats)
    Q_post_modified = deepcopy(Q_post_duplication_normalized)
    Q_pre_modified = deepcopy(Q_pre_duplication)
    nonzerolist = Q_post_duplication_normalized.nonzero()
    nonzeros = [(nonzerolist[0][i],nonzerolist[1][i]) for i in range(len(nonzerolist[0]))]
    for i in range(0,Q_post_duplication_normalized.shape[0]):
        for j in range(0,Q_post_duplication_normalized.shape[1]):
            if not (i,j) in nonzeros:
                Q_post_modified[i,j] += self.err
                
    nonzerolist = Q_pre_duplication.nonzero()
    nonzeros = [(nonzerolist[0][i],nonzerolist[1][i]) for i in range(len(nonzerolist[0]))]
    for i in range(0,Q_pre_duplication.shape[0]):
        for j in range(0,Q_pre_duplication.shape[1]):
            if not (i,j) in nonzeros:
                Q_pre_modified[i,j] += self.err

    # Q modified is just created to avoid 0/0 case        
    C_post = np.divide(Tau_matrix, Q_post_modified) #Coeff Matrix for Geneconv events
    C_post_none_conv = np.divide((Q_post_duplication_normalized - Tau_matrix), Q_post_modified) #Coeff Matrix for none Geneconv events
    np.fill_diagonal(C_post_none_conv,0.0)
    C_pre = np.zeros(Tau_matrix.shape)
    C_pre_none_conv = np.divide(Q_pre_duplication, Q_pre_modified)
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
            edge_to_blen[edge] = self.edge_to_blen[edge]/normalizing_factor
        else:
            edge_to_Q[edge] = Q_pre_duplication
            edge_to_C[edge] = C_pre
            edge_to_C_none[edge] = C_pre_none_conv
            edge_to_blen[edge] = self.edge_to_blen[edge]

    casenum = self.getcasenum(self.modelnum)
    if casenum < 6:
        nt_pairs, pair_to_state = self.get_state_space(2)
        state_pairs = nt_pairs
    else:
        codon_pairs,pair_to_state = self.get_codon_state_space(2)
        state_pairs = codon_pairs

    geneconv_event = 0.0
    none_geneconv_event =0.0
    num_samples = 1000000
    r_list = []
    
    for i in range(num_samples):
        r = test_simulation(Q_pre_duplication,np.zeros(Tau_matrix.shape),d,state_pairs,pair_to_state,edge_to_blen[('N0','Tamarin')])
        geneconv_event += r[2]
        none_geneconv_event += r[3]
        r_list.append(r)

    print geneconv_event/num_samples, none_geneconv_event/num_samples
    Q,d= test.Q_norm(test.modelnum,test.para)

