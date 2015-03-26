from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse

def main(args):
    model = args.model
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
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
        method = 'differential_evolution'
    result = test.get_mle(display = False, em_iterations = 1, method = method)
    test.get_ExpectedNumGeneconv()
    test.get_ExpectedHetDwellTime()
    test.get_individual_summary(summary_path = summary_path)
    test.save_to_file(path = path)

    
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
    
    main(parser.parse_args())    
