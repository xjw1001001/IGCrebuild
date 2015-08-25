from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse

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


    model = 'HKY'
    paralog = ['TLR5a', 'TLR5b']
    alignment_file = './ZebraFish/' + '_'.join(paralog) + '_input.fasta'
    switched_alignment_file = './ZebraFish/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './ZebraFish/FishTree.newick'
    path = './ZebraFish/'
    summary_path = './ZebraFish/'
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
##    result = test.get_mle(display = False, em_iterations = 0, method = method)
##    test.get_ExpectedNumGeneconv()
##    test.get_ExpectedHetDwellTime()
##    test.get_individual_summary(summary_path = summary_path)
##    test.save_to_file(path = path)
