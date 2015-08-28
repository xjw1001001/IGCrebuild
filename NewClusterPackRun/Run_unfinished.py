from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse
import numpy as np
from scipy.linalg import expm
from CodonGeneconFunc import *
import matplotlib.pyplot as plt

def repack_Geneconv_mat(test):
    if test.Model == 'HKY':
        state_size = 4
        mat_basic = test.get_HKYBasic()
        mat_geneconv = np.zeros((state_size**2, state_size**2))
        mat_rcr = test.get_HKYGeneconv()[1]  # rcr stands for row column rate
    elif test.Model == 'MG94':
        state_size = 61
        mat_basic = test.get_MG94Basic()
        mat_geneconv = np.zeros((state_size ** 2, state_size ** 2))
        mat_rcr = test.get_MG94Geneconv()[1]  # rcr stands for row column rate
    row = [ state_size * i[0] + i[1] for i in mat_rcr['row']]
    col = [ state_size * i[0] + i[1] for i in mat_rcr['col']]
    mat_geneconv[row, col] = mat_rcr['rate']

    # now add in diagonal entries
    mat_basic = mat_basic - np.diag(mat_basic.sum(axis = 1))
    mat_geneconv = mat_geneconv - np.diag(mat_geneconv.sum(axis = 1))
    return mat_basic, mat_geneconv
     
def plot_pdiff(test,basic_t, p = None):
    if p == None:
        p = test.prior_distribution
    mat_basic, mat_geneconv = repack_Geneconv_mat(test)
    #basic_t = np.arange(0.0, 1.0, 0.1)
    #geneconv_t = basic_t / 2.0
    basic_psame = [np.matrix(test.prior_distribution) * np.power(expm(mat_basic * t), 2) for t in basic_t]
    if test.Model == 'HKY':
        dup_prior = np.zeros((16))
        identical_state = np.zeros((16))
        for i in range(4):
            dup_prior[i * 4 + i] = test.prior_distribution[i]
            identical_state[i * 4 + i] = 1.0
    elif test.Model == 'MG94':
        dup_prior = np.zeros((61**2))
        identical_state = np.zeros((61**2))
        for i in range(61):
            dup_prior[i * 61 + i] = test.prior_distribution[i]
            identical_state[i * 61 + i] = 1.0

    geneconv_psame = [np.matrix(dup_prior) * np.matrix(expm(mat_geneconv * t)) for t in basic_t]

    basic_pdiff = [1 - t.sum() for t in basic_psame]
    geneconv_pdiff = [1.0 - (identical_state * t.T)[0,0] for t in geneconv_psame]
    
    return basic_pdiff, geneconv_pdiff

def main(args):
    model = args.model
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    switched_alignment_file = './TestTau/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_switched.fasta'
    if args.switch:
        newicktree = './TestTau/YeastTestTree.newick'
        path = './TestTau/'
        summary_path = './TestTau/'
    else:
        newicktree = '../PairsAlignemt/YeastTree.newick'
        path = './NewPackageNewRun/'
        summary_path = './NewPackageNewRun/'
    omega_guess = 0.1    

    print 'Now calculate MLE for pair', paralog
    if args.force:
        if model == 'MG94':
            Force = {5:0.0}
    else:
        Force = None
    if args.gBGC:
        if args.dir:
            test = gBGCDirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
        else:
            test = gBGCCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)       
    else:
        if args.dir:
            test = DirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
        else:
            test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)

    if model == 'MG94':
        method = 'BFGS'
    else:
        method = 'BFGS'
    result = test.get_mle(display = False, em_iterations = 0, method = method)
    test.get_ExpectedNumGeneconv()
    test.get_ExpectedHetDwellTime()
    test.get_individual_summary(summary_path = summary_path)
    test.save_to_file(path = path)


    if args.switch:
        if args.gBGC:
            if args.dir:
                test2 = gBGCDirGeneconv( newicktree, switched_alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
            else:
                test2 = gBGCCodonGeneconv( newicktree, switched_alignment_file, paralog, Model = model, Force = Force, clock = args.clock)       
        else:
            if args.dir:
                test2 = DirGeneconv( newicktree, switched_alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
            else:
                test2 = ReCodonGeneconv( newicktree, switched_alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
        result = test2.get_mle(display = False, em_iterations = 0, method = method)
        test2.get_ExpectedNumGeneconv()
        test2.get_ExpectedHetDwellTime()
        test2.get_individual_summary(summary_path = summary_path + 'switched_')
        test2.save_to_file(path = path + 'switched_')

    
if __name__ == '__main__':
#    parser = argparse.ArgumentParser()
#    parser.add_argument('--model', required = True, help = 'Substitution Model')
#    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
#    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
#    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
#    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
#    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
#    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
#    parser.add_argument('--dir', dest = 'dir', action = 'store_true', help = 'dir control')
#    parser.add_argument('--no-dir', dest = 'dir', action = 'store_false', help = 'dir control')
#    parser.add_argument('--gBGC', dest = 'gBGC', action = 'store_true', help = 'gBGC control')
#    parser.add_argument('--no-gBGC', dest = 'gBGC', action = 'store_false', help = 'gBGC control')
#    parser.add_argument('--switch', dest = 'switch', action = 'store_true', help = 'switch test control')
#    parser.add_argument('--no-switch', dest = 'switch', action = 'store_false', help = 'switch test control')
#    
#    
#    main(parser.parse_args())


##    model = 'HKY'
##    paralog = ['TLR5a', 'TLR5b']
##    alignment_file = './ZebraFish/' + '_'.join(paralog) + '_input.fasta'
##    switched_alignment_file = './ZebraFish/' + '_'.join(paralog) + '_input.fasta'
##    newicktree = './ZebraFish/FishTree.newick'
##    path = './ZebraFish/'
##    summary_path = './ZebraFish/'
##    omega_guess = 0.1
##    force = False
##    gBGC = False
##    Dir = False
##    clock = False
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    #pairs = [pairs[0]]
    
    for paralog in pairs:    
        model = 'HKY'
        #paralog1 = 'YLR406C'
        #paralog2 = 'YDL075W'
        #paralog = [paralog1, paralog2]
        alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        newicktree = './TestTau/YeastTestTree.newick'
        omega_guess = 0.1
        force = False
        gBGC = False
        Dir = False
        clock = False


        print 'Now calculate MLE for pair', paralog
        if force:
            if model == 'MG94':
                Force = {5:0.0}
            else:
                Force = {4:0.0}
        else:
            Force = None
        if gBGC:
            if Dir:
                test = gBGCDirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
            else:
                test = gBGCCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)       
        else:
            if Dir:
                test = DirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
            else:
                test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)

        if model == 'MG94':
            method = 'BFGS'
        else:
            method = 'BFGS'

        test.get_mle(False, True, 0, "BFGS")
        basic_t = np.arange(0.0, 10.0, 0.001)
        basic_pdiff, geneconv_pdiff = plot_pdiff(test, basic_t)
        plt.plot(basic_t, basic_pdiff, 'r-', basic_t, geneconv_pdiff, 'b-')
        plt.ylabel('P_diff')
        plt.title('_'.join(paralog) + ' '+ test.Model + ' long')
        plt.savefig('_'.join(paralog) + ' ' + test.Model + ' long.jpg', bbox_inches='tight')
        plt.close()

        basic_t = np.arange(0.0, 0.2, 0.0001)
        basic_pdiff, geneconv_pdiff = plot_pdiff(test, basic_t)
        plt.plot(basic_t, basic_pdiff, 'r-', basic_t, geneconv_pdiff, 'b-')
        plt.title('_'.join(paralog) + ' '+ test.Model + ' short')
        plt.ylabel('P_diff')
        plt.savefig('_'.join(paralog) + ' ' + test.Model + ' short.jpg', bbox_inches='tight')
        plt.close()
        
    ##    result = test.get_mle(display = False, em_iterations = 0, method = method)
    ##    test.get_ExpectedNumGeneconv()
    ##    test.get_ExpectedHetDwellTime()
    ##    test.get_individual_summary(summary_path = summary_path)
    ##    test.save_to_file(path = path)
