import numpy as np
import os
from GenerateIndividualDirSummary import get_individual_summary

def summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if force:
        prefix = summary_path + 'Force_Dir_' + model + '_'
    else:
        prefix = summary_path + 'Dir_' + model + '_'

    if clock:
        suffix = '_clock_summary.txt'
    else:
        suffix = '_nonclock_summary.txt'

    label = ''
    for pair in pairs:
        summary = prefix + '_'.join(pair) + suffix
        if os.path.isfile(summary):
            res = np.loadtxt(open(summary, 'r'))
            if len(np.atleast_1d(res)) > 1:
                summary_mat.append(res.tolist())
                label = open(summary, 'r').readlines()[-1][2 : -1]
                finished_list.append(pair)
            else:
                unfinished_list.append(pair)
                      
        else:
            unfinished_list.append(pair)
    
    with open(unfinished_list_file, 'w+') as g:
        for pair in unfinished_list:
            g.write('_'.join(pair) + '\n')

    t = np.matrix(summary_mat)
    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    footer = label  # row labels
    np.savetxt(open(summary_file, 'w+'), t.T, delimiter = ' ', header = header, footer = footer)

def get_mass_summary(pairs, pair_path, model, summary_path, clock = True, force = False):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if force:
        prefix = pair_path + 'Force_Dir_' + model + '_'
    else:
        prefix = pair_path + 'Dir_' + model + '_'

    if clock:
        suffix = '_clock.p'
    else:
        suffix = '_nonclock.p'

    label = []
    for pair in pairs:
        p_file = prefix + '_'.join(pair) + suffix
        if os.path.isfile(p_file):
            get_individual_summary(pair, pair_path, model, summary_path, clock, force)


if __name__ == '__main__':
    summary_path = '/Users/xji3/MixedFromCluster/NewPackageNewRun/'
#    summary_path = '/Users/xji3/FromCluster05082015/NewPackageNewRun04102015/'
    summary_path_list=[
        '/Users/xji3/FromCluster03162015/NewPackageNewRun/',
        '/Users/xji3/FromCluster03112015/NewPackageNewRun/',
        '/Users/xji3/FromCluster05082015/NewPackageNewRun04102015/',
        "/Users/xji3/FromCluster03232015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03232015_2/NewPackageNewRun/",
        "/Users/xji3/FromCluster03232015_3/NewPackageNewRun/",
        "/Users/xji3/FromCluster03232015_4/NewPackageNewRun/",
        "/Users/xji3/FromCluster03212015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster05082015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster04172015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03232015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03192015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03182015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03172015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03162015/NewPackageNewRun/",
        "/Users/xji3/FromCluster03112015/NewPackageNewRun/"
##        "/Users/xji3/FromCluster03102015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster03092015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster02192015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster02112015/NewPackageNewRun/",
##        "/Users/xji3/02062015FromCluster/NewPackageNewRun/",
##        "/Users/xji3/FromCluster01272015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster01262015/NewPackageNewRun/",
##        "/Users/xji3/FromCluster01162015/NewPackageNewRun/",
##        "/Users/xji3/0112FromCluster/NewPackageNewRun/",
##        "/Users/xji3/1229FromCluster/NewPackageNewRun/"
        ]
    
    pairs = []
    all_pairs = '../All_Pairs.txt'
    jeff_pairs = './Jeff_pairs_list.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])
    summary_path_list = ['/Users/xji3/MixedFromCluster/NewPackageNewRun/']
    summary_path_list = ["/Users/xji3/FromCluster06082015/TestTau/", '/Users/xji3/FromCluster06082015/TestTau/switched_']
####################################################################################################################################################
##
    for summary_path in summary_path_list:
        # MG94 clock model
        model = 'MG94'
        summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_MG94_clock_summary.txt',
                         unfinished_list_file = summary_path + 'Dir_MG94_clock_unfinished.txt',
                         clock = True, force = False)

        # MG94 nonclock model 
        summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_MG94_nonclock_summary.txt',
                         unfinished_list_file = summary_path + 'Dir_MG94_nonclock_unfinished.txt',
                         clock = False, force = False)

    ##    # MG94 clock model (force tau = 0)
    ##    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_Force_MG94_clock_summary.txt',
    ##                     unfinished_list_file = summary_path + 'Dir_Force_MG94_clock_unfinished.txt',
    ##                     clock = True, force = True)
    ##    # MG94 nonclock model (force tau = 0)
    ##    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_Force_MG94_nonclock_summary.txt',
    ##                     unfinished_list_file = summary_path + 'Dir_Force_MG94_nonclock_unfinished.txt',
    ##                     clock = False, force = True)
    ##                

    ####################################################################################################################################################

        model = 'HKY'
    ##    # HKY clock model
    ##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
    ##                     clock = True, force = False)
    ##
    ##    # HKY nonclock model 
    ##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
    ##                     clock = False, force = False)

        # HKY clock model
        summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_HKY_clock_summary.txt',
                         unfinished_list_file = summary_path + 'Dir_HKY_clock_unfinished.txt',
                         clock = True, force = False)

        # HKY nonclock model 
        summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_HKY_nonclock_summary.txt',
                         unfinished_list_file = summary_path + 'Dir_HKY_nonclock_unfinished.txt',
                         clock = False, force = False)
