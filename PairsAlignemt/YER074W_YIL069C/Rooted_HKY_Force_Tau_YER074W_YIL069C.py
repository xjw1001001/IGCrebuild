'''
Gene Pair YER074W YIL069C
'''
import sys
sys.path.append('../')
from CodonBased2Repeats import *

if __name__=='__main__':
    numLeaf = 7
    blen = np.ones([2*numLeaf-2])*2

    tree_newick = '../YeastTree.newick'
    dataloc = './YER074W_YIL069C_input.fasta'
    paralog_pair = ['YER074W','YIL069C']

    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,cdmodel = False, removegaps = True)

    Model = 'HKY'
    guess = np.log([0.49063051,   0.52677502,   0.37278418,   0.9])
    #guess = np.append(guess,1.5/2.0)

    leaves = set(v for v, degree in test.treetopo.degree().items() if degree == 1)
    guess_w_blen = [-2.3]*(len(leaves))
    guess_w_blen.extend(guess)

    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()
    
    r1 = test.estimate(args, Model, guess_w_blen,output = './RootedTest_HKY'+'_'.join(test.paralog)+'_Force_Tau.txt', est_blen = True,force_tau = True, clock = True,print_result = True)

    test.save_to_file('./YER074W_YIL069C_R0_HKY_Class.p')

