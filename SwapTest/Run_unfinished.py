from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse
from Pdiff import *


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
        elif model == 'HKY':
            Force = {4:0.0}
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
    if model == 'MG94':
        test.get_individual_summary(summary_path = summary_path)


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


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required = True, help = 'Substitution Model')
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
    parser.add_argument('--dir', dest = 'dir', action = 'store_true', help = 'dir control')
    parser.add_argument('--no-dir', dest = 'dir', action = 'store_false', help = 'dir control')
    parser.add_argument('--gBGC', dest = 'gBGC', action = 'store_true', help = 'gBGC control')
    parser.add_argument('--no-gBGC', dest = 'gBGC', action = 'store_false', help = 'gBGC control')
    parser.add_argument('--switch', dest = 'switch', action = 'store_true', help = 'switch test control')
    parser.add_argument('--no-switch', dest = 'switch', action = 'store_false', help = 'switch test control')
    
    
    main(parser.parse_args())


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


##    pairs = []
##    all_pairs = './Filtered_pairs.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##    #pairs = [pairs[0]]
##
####
####    pairs = [['EDN', 'ECP']]
####    
##    for paralog in pairs:    
##        model = 'MG94'
##        #model = 'HKY'
##        #paralog1 = 'YLR406C'
##        #paralog2 = 'YDL075W'
##        #paralog = [paralog1, paralog2]
##        alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
##        #alignment_file = '../data/EDN_ECP_Cleaned.fasta'
##        newicktree = './TestTau/YeastTestTree.newick'
##        newicktree = '../PairsAlignemt/YeastTree.newick'
##        #newicktree = '../data/input_tree.newick'
##        
##        omega_guess = 0.1
##        force = False
##        gBGC = False
##        Dir = False
##        clock = False
##        summary_path = '/Users/xji3/MixedFromCluster/NewPackageNewRun/'
##        #summary_path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/OldResults02112015/'
##
##
##        print 'Now calculate MLE for pair', paralog
##        if force:
##            if model == 'MG94':
##                Force = {5:0.0}
##            else:
##                Force = {4:0.0}
##        else:
##            Force = None
##        if gBGC:
##            if Dir:
##                test = gBGCDirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
##            else:
##                test = gBGCCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)       
##        else:
##            if Dir:
##                test = DirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
##            else:
##                test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = clock)
##
##        if model == 'MG94':
##            method = 'BFGS'
##        else:
##            method = 'BFGS'
##        x = read_txt(summary_path, paralog, model, force, clock, Dir, gBGC)#, Spe = 'Primate')
##        test.update_by_x(x)
##
##        new_summary_path = './NewPackageNewRun/'
##        test.get_individual_summary(new_summary_path)
        #test.get_mle(True, True, 0, method)

##        #### Get plots for specific range
##        max_point = 10.0
##        number_dots = 100
##        step = max_point / number_dots 
##        basic_t = np.arange(0.0, max_point, step)
##        pdiff = np.loadtxt(open('/Users/xji3/plotPdiff09032015/' + '_'.join(paralog) + ' ' + test.Model + ' max ' + str(max_point) + ' data.txt'))
##        basic_pdiff_short, geneconv_pdiff_short = pdiff[:number_dots], pdiff[number_dots:]
####        basic_pdiff_short, geneconv_pdiff_short, mut_odds_short, geneconv_odds_short = plot_pdiff(test, basic_t)
##        if abs(basic_pdiff_short[0]) < 1e-10:
##            basic_pdiff_short[0] = 0.0
##        if abs(geneconv_pdiff_short[0]) < 1e-10:
##            geneconv_pdiff_short[0] = 0.0
####        np.savetxt(open('_'.join(paralog) + ' ' + test.Model + ' max ' + str(max_point) + ' data.txt', 'w+'),np.concatenate((basic_pdiff_short, geneconv_pdiff_short), axis = 0))
##        plt.plot(2 * basic_t, basic_pdiff_short, 'r-', label = 'Basic Model')
##        plt.plot(2 * basic_t, geneconv_pdiff_short, 'b-', label = 'IGC Model')
##        plt.legend(bbox_to_anchor = (1,1), loc = 2, borderaxespad = 0.)
##        plt.ylabel('P_diff')
##        plt.title('_'.join(paralog) + ' '+ test.Model + ' max ' + str(max_point))
##        plt.savefig('_'.join(paralog) + ' ' + test.Model + ' max ' + str(max_point) + '.pdf', bbox_inches='tight')
##        plt.close()

##        Mut = np.loadtxt(open('/Users/xji3/plotPdiff09032015/' + '_'.join(paralog) + ' ' + test.Model + ' Mut IGC max ' + str(max_point) + ' data.txt'))
##        mut_odds_short, geneconv_odds_short = Mut[:number_dots], Mut[number_dots:]
##        if abs(mut_odds_short[0]) < 1e-10:
##            mut_odds_short[0] = 0.0
##        if abs(geneconv_odds_short[0]) < 1e-10:
##            geneconv_odds_short[0] = 0.0
####        np.savetxt(open('_'.join(paralog) + ' ' + test.Model + ' Mut IGC max ' + str(max_point) + ' data.txt', 'w+'),np.concatenate((mut_odds_short, geneconv_odds_short), axis = 0))
####        
##        plt.plot(basic_t, mut_odds_short, 'r-', label = 'Mutation Odds')
##        plt.plot(basic_t, geneconv_odds_short, 'b-', label = 'IGC Odds')
##        plt.legend(bbox_to_anchor = (1,1), loc = 2, borderaxespad = 0.)
##        plt.ylabel('Odds')
##        plt.title('_'.join(paralog) + ' '+ test.Model + ' Mut IGC max ' + str(max_point))
##        plt.savefig('_'.join(paralog) + ' ' + test.Model + ' Mut IGC max ' + str(max_point) + '.jpg', bbox_inches='tight')
##        plt.close()


        

    ##    result = test.get_mle(display = False, em_iterations = 0, method = method)
    ##    test.get_ExpectedNumGeneconv()
    ##    test.get_ExpectedHetDwellTime()
    ##    test.get_individual_summary(summary_path = summary_path)
    ##    test.save_to_file(path = path)
