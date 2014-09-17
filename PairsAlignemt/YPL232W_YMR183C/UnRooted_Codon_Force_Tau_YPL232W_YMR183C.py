'''
Gene Pair YPL232W YMR183C
'''
import sys
sys.path.append('../')
from CodonBased2Repeats import *

if __name__=='__main__':
    numLeaf = 7
    blen = np.ones([2*numLeaf-2])*2

    tree_newick = '../YeastTree.newick'
    dataloc = './YPL232W_YMR183C_input.fasta'
    paralog_pair = ['YPL232W','YMR183C']

    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,cdmodel = True, removegaps = True)

    Model = 'Codon_corrected'
    guess = np.log([0.49063051,   0.52677502,   0.37278418,   8.45263599,  17.77868809])
    #guess = np.append(guess,0.191996871131)

    guess_w_blen = np.log([0.06767592097913287, 0.008039951903822255, 0.06452861394808071, 0.1274028147853794, 0.2635103777495246, 0.00802702317557723, 0.053738461527140756, 0.09825885747458973, 0.07306103166567632, 0.09458992499865115, 0.1881915684864276, 0.02482047743757791])
    guess_w_blen = np.append(guess_w_blen, guess)

    parser = argparse.ArgumentParser()
    parser.add_argument('--n', default=10, type=int,
            help='number of sites in the sampled alignment')
    parser.add_argument('--phi', default=2.0, type=float,
            help='strength of the gene conversion effect')
    args = parser.parse_args()
    
    r1 = test.estimate(args, Model, guess_w_blen,output = './UnRootedTest_Codon_corrected'+'_'.join(test.paralog)+'_Force_Tau.txt',force_tau = True, est_blen = True, clock = False,print_result = True)


