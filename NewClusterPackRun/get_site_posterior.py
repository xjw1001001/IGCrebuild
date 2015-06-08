from Rewrite_CodonGeneconv import ReCodonGeneconv
from DirGeneconv import DirGeneconv
from gBGCDirGeneconv import gBGCDirGeneconv
from gBGCCodonGeneconv import gBGCCodonGeneconv
from CodonGeneconFunc import *
import argparse
import numpy as np
import scipy
import pickle
from copy import deepcopy
from functools import partial

# jsconctmctree commit number f510668605e43c2596eec9fe6f99fb639054112d


def readSummaryFile(summary_file):
    return np.loadtxt(summary_file)

def get_site_ExpectedNumGeneconv(test):
    if test.GeneconvTransRed is None:
        test.GeneconvTransRed = test.get_directionalNumGeneconvRed()

    scene_ll = test.get_scene()
    requests = [{'property' : 'DDNTRAN', 'transition_reduction' : i} for i in test.GeneconvTransRed]
    assert(len(requests) == 2)
    j_in = {
        'scene' : scene_ll,
        'requests' : requests
        }        
    j_out = jsonctmctree.interface.process_json_in(j_in)

    status = j_out['status']
    ExpectedGeneconv_mat12 = np.matrix(j_out['responses'][0]).T
    ExpectedGeneconv_mat21 = np.matrix(j_out['responses'][1]).T
    #ExpectedGeneconv = {test.edge_list[i] : ExpectedGeneconv_mat[i,:].tolist()[0] for i in range(len(test.edge_list))}
    return ExpectedGeneconv_mat12, ExpectedGeneconv_mat21

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
        
    test.x = transformed_x
    test.update_by_x()
        
    #b = pickle.load(open(prefix + '_'.join(paralog) + suffix_pickle, 'rb'))
    
    test._loglikelihood2()
    print test.ll, gBGC#, b['ll']#, transformed_x

    ExpectedGeneconv_mat12, ExpectedGeneconv_mat21 = get_site_ExpectedNumGeneconv(test)
    footer = ' '.join(['_'.join(edge) + '_12' for edge in test.edge_list]) + ' '.join(['_'.join(edge) + '_21' for edge in test.edge_list])
    np.savetxt(open(prefix.replace(summary_path, path) + '_'.join(paralog) + suffix.replace('summary', 'SiteExpectedGeneconv'), 'w+'),
               np.concatenate((ExpectedGeneconv_mat12, ExpectedGeneconv_mat21)), delimiter = ' ', footer = footer)    
    

if __name__ == '__main__':
   
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--model', required = True, help = 'Substitution Model')
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--force', dest = 'force', action = 'store_true', help = 'Tau parameter control')
##    parser.add_argument('--no-force', dest = 'force', action = 'store_false', help = 'Tau parameter control')
##    parser.add_argument('--clock', dest = 'clock', action = 'store_true', help = 'clock control')
##    parser.add_argument('--no-clock', dest = 'clock', action = 'store_false', help = 'clock control')
##    parser.add_argument('--dir', dest = 'dir', action = 'store_true', help = 'dir control')
##    parser.add_argument('--no-dir', dest = 'dir', action = 'store_false', help = 'dir control')
##    parser.add_argument('--gBGC', dest = 'gBGC', action = 'store_true', help = 'gBGC control')
##    parser.add_argument('--no-gBGC', dest = 'gBGC', action = 'store_false', help = 'gBGC control')
##    
##    main(parser.parse_args())    

    paralog = ['YDR418W', 'YEL054C']
    model = 'HKY'
    force = False
    clock = False
    gBGC = False
    Dir = False
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
        
    test.x = transformed_x
    test.update_by_x()
        
    #b = pickle.load(open(prefix + '_'.join(paralog) + suffix_pickle, 'rb'))
    
    test._loglikelihood2()
    print test.ll, gBGC#, b['ll']#, transformed_x

    ExpectedGeneconv_mat12, ExpectedGeneconv_mat21 = get_site_ExpectedNumGeneconv(test)
    footer = ' '.join(['_'.join(edge) + '_12' for edge in test.edge_list]) + ' '.join(['_'.join(edge) + '_21' for edge in test.edge_list])
    np.savetxt(open(prefix.replace(summary_path, path) + '_'.join(paralog) + suffix.replace('summary', 'SiteExpectedGeneconv'), 'w+'),
               np.concatenate((ExpectedGeneconv_mat12, ExpectedGeneconv_mat21)), delimiter = ' ', footer = footer)    
     
    



