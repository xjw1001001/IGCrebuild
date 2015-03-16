import os.path
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import pickle
import numpy as np
import argparse

def get_summary(p_file, directional = False, output_label = False):
    res = pickle.load(open(p_file, 'r'))
    if directional:
        t = gBGCDirGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'], clock = res['clock'], Force = res['Force'])
    else:
        t = gBGCCodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'], clock = res['clock'], Force = res['Force'])
        
    if res['clock']:
        t.update_by_x_clock(res['x_clock'])
    else:
        t.update_by_x(res['x'])

   
    out = [t.nsites, res['ll']]
    out.extend(t.pi)
    if directional:
        if t.Model == 'HKY': # HKY model doesn't have omega parameter
            out.extend([t.kappa])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau12', 'tau21', 'gamma']
        elif t.Model == 'MG94':
            out.extend([t.kappa, t.omega])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau12', 'tau21', 'gamma']
        out.extend(t.tau)
        out.extend([t.gamma])
    else:
        if t.Model == 'HKY': # HKY model doesn't have omega parameter
            out.extend([t.kappa, t.tau, t.gamma])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau', 'gamma']
        elif t.Model == 'MG94':
            out.extend([t.kappa, t.omega, t.tau, t.gamma])
            label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau', 'gamma']
        
    k = len(label)  # record the length of non-blen parameters

    label.extend(res['edge_list'])

    out.extend([t.edge_to_blen[label[j]] for j in range(k, len(label))])

    if res.has_key('ExpectedGeneconv'):
        t.ExpectedGeneconv = res['ExpectedGeneconv']
    if res.has_key('ExpectedDwellTime'):
        t.ExpectedDwellTime = res['ExpectedDwellTime']

    if not t.ExpectedGeneconv:
        t.get_ExpectedNumGeneconv()

    if not t.ExpectedDwellTime:
        t.get_ExpectedHetDwellTime()

    label.extend([ (a, b, 'tau') for (a, b) in res['edge_list']])
    if directional:
        out.extend([sum(t.ExpectedGeneconv[i]) / (t.edge_to_blen[i] * t.ExpectedDwellTime[i]) if t.ExpectedDwellTime[i] != 0 else 0 for i in res['edge_list']])
    else:
        out.extend([t.ExpectedGeneconv[i] / (t.edge_to_blen[i] * t.ExpectedDwellTime[i]) if t.ExpectedDwellTime[i] != 0 else 0 for i in res['edge_list']])


    # Now add directional # of geneconv events
    ExpectedDirectionalNumGeneconv = t._ExpectedDirectionalNumGeneconv()
    label.extend([ (a, b, '1->2') for (a, b) in res['edge_list']])
    out.extend([ExpectedDirectionalNumGeneconv[i][0] for i in res['edge_list']])
    label.extend([ (a, b, '2->1') for (a, b) in res['edge_list']])
    out.extend([ExpectedDirectionalNumGeneconv[i][1] for i in res['edge_list']])


    for i in range(k, len(label)):
        label[i] = '(' + ','.join(label[i]) + ')'

    if output_label:
        return out, label
    else:
        return out

def get_individual_summary(pair, pair_path, model, summary_path, clock, force, directional):
    if directional:
        if force:
            prefix = pair_path + 'Force_gBGC_Dir_' + model + '_'
            prefix_summary = summary_path + 'Force_gBGC_Dir_' + model + '_'
        else:
            prefix = pair_path + 'gBGC_Dir_' + model + '_'
            prefix_summary = summary_path + 'gBGC_Dir_' + model + '_'
    else:        
        if force:
            prefix = pair_path + 'Force_gBGC_' + model + '_'
            prefix_summary = summary_path + 'Force_gBGC_' + model + '_'
        else:
            prefix = pair_path + 'gBGC_' + model + '_'
            prefix_summary = summary_path + 'gBGC_' + model + '_'

    if clock:
        suffix = '_clock.p'
        suffix_summary = '_clock_summary.txt'
    else:
        suffix = '_nonclock.p'
        suffix_summary = '_nonclock_summary.txt'    

    p_file = prefix + '_'.join(pair) + suffix
    summary_file = prefix_summary + '_'.join(pair) + suffix_summary
    if os.path.isfile(p_file):
        res = get_summary(p_file, directional, True)
        summary = np.matrix(res[0])
        label = res[1]
    else:
        summary = np.matrix(-1.234)
        label = ['unfinished']
        
    footer = ' '.join(label)  # row labels
    np.savetxt(open(summary_file, 'w+'), summary.T, delimiter = ' ', footer = footer)

def main(args):
    pair = [args.paralog1, args.paralog2]
    clock = args.clock
    force = args.force
    directional = args.dir
    summary_path = args.sump
    model = args.model
    pair_path = args.pairp

    get_individual_summary(pair, pair_path, model, summary_path, clock, force, directional)
    
if __name__ == '__main__':
    
##    p_file = '/Users/xji3/FromCluster03162015/NewClusterPackRun/NewPackageNewRun/gBGC_HKY_YAL056W_YOR371C_clock.p'
##    p_file = '/Users/xji3/FromCluster03162015/NewClusterPackRun/NewPackageNewRun/gBGC_Dir_HKY_YAL056W_YOR371C_clock.p'
##    print get_summary(p_file, True)
##    pair = ['YAL056W', 'YOR371C']
##    get_individual_summary(pair, '/Users/xji3/FromCluster03162015/NewClusterPackRun/NewPackageNewRun/', 'HKY', './', False, None, False)
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
    parser.add_argument('--model', required = True, help = 'model control')
    parser.add_argument('--sump', required = True, help = 'summary file path')
    parser.add_argument('--pairp', required = True, help = 'pair file path')
    parser.add_argument('--dir', dest = 'dir', action = 'store_true', help = 'dir control')
    parser.add_argument('--no-dir', dest = 'dir', action = 'store_false', help = 'dir control')
    main(parser.parse_args())
