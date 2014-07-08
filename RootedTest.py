from CodonBased2Repeats import *

if __name__ == '__main__':
    #sys.stdout = open('RootedTestPrintoScreen.txt','w+')
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
##
    leaves = set(v for v, degree in test.treetopo.degree().items() if degree == 1) 
    guess_w_blen = [-1.0]*(len(leaves))
    #guess_w_blen[0:6] = [0.0,-1.0,-2.0,-3.0,-4.0,-5.0]
    guess_w_blen.extend(guess)
    r1 = test.estimate(args, Model, guess_w_blen,output = './RootedTest_HKY.txt', est_blen = True, clock = True,print_result = False)

        

    


    
    numLeaf_codon = 5
    blen_codon = np.ones([2*numLeaf_codon-2])*2
    dataloc_codon = './data/codon_alignment_nucleotide_format.fasta'
    tree_newick = './data/input_tree.newick'
    Model = 'Codon_corrected'
    test4  = Codon2RepeatsPhy(numLeaf_codon,blen_codon,tree_newick,dataloc_codon,cdmodel= True,removegaps = True)
    leaves = set(v for v, degree in test4.treetopo.degree().items() if degree == 1)
    guess_w_blen = [-1.0]*(len(leaves))
    guess = np.log([0.49628073,  0.58854499 , 0.48857723,  2.10157473,  1.16726338])
    guess_w_blen = np.append(guess_w_blen, guess)
    r5 = test4.estimate(args, 7, guess_w_blen,output = './RootedTest_Codon_corrected.txt', est_blen = True,clock = True,force_tau = True, print_result = True)
    
