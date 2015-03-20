import numpy as np
import os
from GenerateIndividualgBGCSummary import get_individual_summary

def get_mass_summary(pairs, pair_path, model, summary_path, clock, force, directional):
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
            get_individual_summary(pair, pair_path, model, summary_path, clock, force, directional)



def summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, directional):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if directional:
        if force:
            prefix = summary_path + 'Force_gBGC_Dir_' + model + '_'
        else:
            prefix = summary_path + 'gBGC_Dir_' + model + '_'
    else:        
        if force:
            prefix = summary_path + 'Force_gBGC_' + model + '_'
        else:
            prefix = summary_path + 'gBGC_' + model + '_'

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
            else:
                unfinished_list.append(pair)
                
            finished_list.append(pair)
        else:
            unfinished_list.append(pair)
    
    with open(unfinished_list_file, 'w+') as g:
        for pair in unfinished_list:
            g.write('_'.join(pair) + '\n')

    t = np.matrix(summary_mat)
    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    footer = label  # row labels
    np.savetxt(open(summary_file, 'w+'), t.T, delimiter = ' ', header = header, footer = footer)

if __name__ == '__main__':
    summary_path = '/Users/xji3/FromCluster03192015/NewPackageNewRun/'
    model = 'MG94'
    pairs = []
    all_pairs = '../All_Pairs.txt'
    jeff_pairs = './Jeff_pairs_list.txt'
    with open(jeff_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])

####################################################################################################################################################

    # MG94 clock model
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_MG94_clock_summary.txt',
                     unfinished_list_file = summary_path + 'gBGC_MG94_clock_unfinished.txt',
                     clock = True, force = False, directional = False)

    # MG94 nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_MG94_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'gBGC_MG94_nonclock_unfinished.txt',
                     clock = False, force = False, directional = False)

    # MG94 directional clock model
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_MG94_clock_summary.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_MG94_clock_unfinished.txt',
                     clock = True, force = False, directional = True)

    # MG94 directional nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_MG94_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_MG94_nonclock_unfinished.txt',
                     clock = False, force = False, directional = True)

####################################################################################################################################################

    model = 'HKY'
##    # HKY clock model
##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
##                     clock = True, force = False, directional = False)
##
##    # MG94 nonclock model 
##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
##                     clock = False, force = False, directional = False)
##
##    # MG94 directional clock model
##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
##                     clock = True, force = False, directional = True)
##
##    # MG94 directional nonclock model 
##    get_mass_summary(pairs, summary_path, model, summary_path = summary_path,
##                     clock = False, force = False, directional = True)

    # HKY clock model
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_HKY_clock_summary.txt',
                     unfinished_list_file = summary_path + 'gBGC_HKY_clock_unfinished.txt',
                     clock = True, force = False, directional = False)

    # HKY nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_HKY_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'gBGC_HKY_nonclock_unfinished.txt',
                     clock = False, force = False, directional = False)

    # HKY directional clock model
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_HKY_clock_summary.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_HKY_clock_unfinished.txt',
                     clock = True, force = False, directional = True)

    # HKY directional nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_HKY_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_HKY_nonclock_unfinished.txt',
                     clock = False, force = False, directional = True)
