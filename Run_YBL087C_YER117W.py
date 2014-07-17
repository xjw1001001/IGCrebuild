'''
Gene Pair YBL087C YER117W
'''
from CodonBased2Repeats import *

if __name__=='__main__':
    numLeaf = 7
    blen = np.ones([2*numLeaf-2])*2

    tree_newick = './data/YeastTree.newick'
    dataloc = './data/Pairs/YBL087C_YER117W/YBL087C_YER117W_input_Aligned.fasta'
    paralog_pair = ['YBL087C','YER117W']

    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,align=False, removegaps = True)

    Model = 'HKY'
    guess = np.log([0.7,0.5,0.4,np.e**2])
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
    
    r1 = test.estimate(args, Model, guess_w_blen,output = './RootedTest_HKY'+'_'.join(test.paralog)+'_Free_Tau.txt', est_blen = True, clock = True,print_result = False)


    test2 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,align=False, removegaps = True)
    guess = np.log([0.7,0.5,0.4,np.e**2])
    guess_w_blen = [-1.0]*(len(leaves))
    guess_w_blen.extend(guess)
    r2 = test2.estimate(args, Model, guess_w_blen,output = './RootedTest_HKY'+'_'.join(test.paralog)+'_Force_Tau.txt', est_blen = True, clock = True,force_tau = True,print_result = False)


    test3 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,align=False, removegaps = True)
    guess_w_blen = np.log([0.27833287642273247, 0.2859705023574136, 0.27718461084282914, 0.3085589061697145, 0.2826254018176859, 0.30855555299681786, 0.2903408485461482, 0.29567177183642623, 0.28490413150490823, 0.2981822976574691, 0.2843326097813424, 0.28780061271138896])
    guess = np.log([0.7,0.5,0.4,np.e**2])
    guess = np.append(guess,14.0)
    guess_w_blen = np.append(guess_w_blen, guess)
    r3 = test3.estimate(args, Model, guess_w_blen,output = './Unrooted_HKY'+'_'.join(test.paralog)+'_Free_Tau.txt',unrooted = False, est_blen = True,clock = False,force_tau = False, print_result = False)

    test4 = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,align=False, removegaps = True)
    guess_w_blen = np.log([0.27833287642273247, 0.2859705023574136, 0.27718461084282914, 0.3085589061697145, 0.2826254018176859, 0.30855555299681786, 0.2903408485461482, 0.29567177183642623, 0.28490413150490823, 0.2981822976574691, 0.2843326097813424, 0.28780061271138896])
    guess = np.log([0.7,0.5,0.4,np.e**2])
    #guess = np.append(guess,14.0)
    guess_w_blen = np.append(guess_w_blen, guess)
    r4 = test4.estimate(args, Model, guess_w_blen,output = './Unrooted_HKY'+'_'.join(test.paralog)+'_Force_Tau.txt',unrooted = False, est_blen = True,clock = False,force_tau = True, print_result = False)


    
