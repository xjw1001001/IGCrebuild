from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse
import numpy as np

def main(args):
    model = args.model
    paralog = [args.paralog1, args.paralog2]
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    path = './NewPackageNewRun/'
    summary_path = './NewPackageNewRun/'
    omega_guess = 0.1    

    Force = None

    txtname = '_'.join(paralog) + '_' + model

    if args.gBGC:
        txtname = txtname + '_gBGC'
        if args.dir:
            test = gBGCDirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
            txtname = txtname + '_dir'
        else:
            test = gBGCCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)       
    else:
        if args.dir:
            txtname = txtname + '_dir'
            test = DirGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)
        else:
            test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = model, Force = Force, clock = args.clock)

    if test.clock:
        txtname = txtname + '_clock'
    else:
        txtname = txtname + '_nonclock'
        
    Force_index = len(test.x_process) + 1

    out_group_blen = np.arange(0.0001, 0.01, 0.0005)
    ll_list = []
    for blen in out_group_blen:
        test.Force = {Force_index:blen}
        test.get_initial_x_process()
        test.update_by_x()
        test.get_mle(False, True, 1, 'BFGS')
        ll_list.append(test.ll)


    np.savetxt(open('./testlikelihood_' + txtname + '.txt', 'w+'), np.matrix([ll_list, out_group_blen]), delimiter = ' ')

    
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
