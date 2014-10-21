'''
Gene Pair YBR024W YBR037C
'''
import sys
sys.path.append('../')
from CodonBased2Repeats import *

if __name__=='__main__':
    numLeaf = 7
    blen = np.ones([2*numLeaf-2])*2

    tree_newick = '../YeastTree.newick'
    dataloc = './YBR024W_YBR037C_input.fasta'
    paralog_pair = ['YBR024W','YBR037C']

    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,cdmodel = False, removegaps = True)

    Model = 'HKY'
    guess = np.log([0.49628073,  0.58854499 , 0.48857723,  2.10157473])
    guess = np.append(guess,1.5/2.0)

    leaves = set(v for v, degree in test.treetopo.degree().items() if degree == 1)
    guess_w_blen = [-1.0]*(len(leaves))
    guess_w_blen.extend(guess)

    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()
    
    r1 = test.estimate(args, Model, guess_w_blen,output = './RootedTest_HKY'+'_'.join(test.paralog)+'_Free_Tau.txt', est_blen = True, clock = True,print_result = True)

    test.save_to_file('./YBR024W_YBR037C_RF_HKY_Class.p')

