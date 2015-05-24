from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
import argparse
import numpy as np
import scipy
import pickle
from copy import deepcopy
from functools import partial

import numdifftools as nd

def readSummaryFile(summary_file):
    return np.loadtxt(summary_file)
    
def get_hessian(test,gBGC = False, epsilon = 1e-5):
    if not test.clock:
        xx_process = deepcopy(test.x_process)
        xx_rates = deepcopy(test.x_rates)
        f = nd.Hessian(partial(get_likelihood_nonclock, test = test, gBGC = gBGC))
        if gBGC:
            xx_process[:-1] = np.exp(xx_process[:-1])
        else:
            xx_process = np.exp(xx_process)
        xx_process[0:3] = np.log(-np.log(xx_process[0:3]))
        xx_rates = np.exp(xx_rates)
        hessian_process = f(xx_process)
        hessian_rates = np.zeros((len(xx_rates), len(test.x)))
        for j in range(len(xx_rates)):
            xx_rates[j] += epsilon / 2.0
            test.update_by_x(np.concatenate((test.x_process, np.log(xx_rates))))
            f2 = -test.loglikelihood_and_gradient()[1] /np.concatenate((np.exp(test.x_process), xx_rates))
            xx_rates[j] -= epsilon
            test.update_by_x(np.concatenate((test.x_process, np.log(xx_rates))))
            f1 = -test.loglikelihood_and_gradient()[1] / np.concatenate((np.exp(test.x_process), xx_rates))
            xx_rates[j] += epsilon / 2.0
            test.update_by_x(np.concatenate((test.x_process, np.log(xx_rates))))
            hessian_rates[j, :] = (f2 - f1) / epsilon

        hessian = np.zeros((len(test.x), len(test.x)))
        hessian[:len(test.x_process), :len(test.x_process)] = hessian_process
        hessian[len(test.x_process):len(test.x), :] = hessian_rates
        hessian[:len(test.x_process), len(test.x_process):] = hessian_rates[:, :len(test.x_process)].T

    return hessian

def get_likelihood_nonclock(x, test, gBGC = False):
    assert(not test.clock)
    #print x
    x[0:3] = -np.exp(x[0:3]) # log(-log) transform of first three elements
    if gBGC:
        x[3:-1] = np.log(x[3:-1])
    else:
        x[3:] = np.log(x[3:])    # no transform of other elements
    test.update_by_x(np.concatenate((x, test.x_rates)))
    return test._loglikelihood2()[0]

        

def get_hessian_BF(test,gBGC = False, epsilon = 1e-5):
    if not test.clock:
        xx = deepcopy(test.x)
        num_rates = len(test.x_rates)
        f = nd.Hessian(partial(get_likelihood_nonclock_BF, test = test, gBGC = gBGC))
        if gBGC:
            xx[:-(1 + num_rates)] = np.exp(xx[:-(1 + num_rates)])
        else:
            xx[:-num_rates] = np.exp(xx[:-num_rates])
        xx[0:3] = np.log(-np.log(xx[0:3]))
        xx[-num_rates:] = np.exp(xx[-num_rates:])
        hessian = f(xx)

    return hessian         

def get_likelihood_nonclock_BF(x, test, gBGC = False):
    assert(not test.clock)
    #print x
    x[0:3] = -np.exp(x[0:3]) # log(-log) transform of first three elements
    num_rates = len(test.x_rates)
    if gBGC:
        x[3:-(1 + num_rates)] = np.log(x[3:-(1 + num_rates)])
    else:
        x[3: -num_rates] = np.log(x[3: -num_rates])    # no transform of other elements
    x[-num_rates :] = np.log(abs(x[-num_rates :]))
    test.update_by_x(x)
    return test._loglikelihood2()[0]

def main(args):
    paralog = [args.paralog1, args.paralog2]
    model = args.model
    force = args.force
    clock = args.clock
    gBGC = args.gBGC
    Dir = args.dir
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = '../PairsAlignemt/YeastTree.newick'
    path = './NewPackageNewRun/'
    summary_path = "../MixedFromCluster/NewPackageNewRun/"


    if force:
        prefix = summary_path + 'Force_' + model + '_'
    else:
        if Dir:
            if gBGC:
                prefix = summary_path + 'gBGC_Dir_' + model + '_'
            else:
                prefix = summary_path + 'Dir_' + model + '_'
        else:
            if gBGC:
                prefix = summary_path + 'gBGC_' + model + '_'
            else:
                prefix = summary_path + model + '_'

    if clock:
        suffix = '_clock_summary.txt'
        suffix_pickle = '_clock.p'
    else:
        suffix = '_nonclock_summary.txt'
        suffix_pickle = '_nonclock.p'
    

    print 'Now calculate Hessian Matrix for pair', paralog
    if force:
        if model == 'MG94':
            Force = {5:0.0}
        elif model == 'HKY':
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


    ind_summary = readSummaryFile(prefix + '_'.join(paralog) + suffix)
    untransformed_x = ind_summary[2:(3 + len(test.x))]
    transformed_x = np.log(np.concatenate(([untransformed_x[0] + untransformed_x[2],
                                        untransformed_x[0] / (untransformed_x[0] + untransformed_x[2]),
                                        untransformed_x[1] / (untransformed_x[1] + untransformed_x[3])],
                                       abs(untransformed_x[4:]))))
    if gBGC:
        transformed_x[len(test.x_process) - 1] = untransformed_x[len(test.x_process)]

    test.update_by_x(transformed_x)
    if clock:
        test.update_x_clock_by_x()
        #test.update_by_x_clock()
        
    b = pickle.load(open(prefix + '_'.join(paralog) + suffix_pickle, 'rb'))
    
    test._loglikelihood2()
    print test.ll, b['ll']
    if clock:
        
        for i in test.edge_list:
            print i, test.edge_to_blen[i]
            
        test.update_x_clock_by_x()
        test.update_by_x_clock()
        for i in test.edge_list:
            print i, test.edge_to_blen[i]

    if model == 'HKY':
        h = get_hessian_BF(test, gBGC, 1e-5)
    else:
        h = get_hessian(test, epsilon = 1e-5, gBGC = gBGC)

    np.savetxt(open(prefix + '_'.join(paralog) + suffix.replace('summary', 'hessian'), 'w+'), h, delimiter = ' ')
    #h_inverse = np.linalg.inv(-h)
    
    

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
   
    
