import numpy as np
import os

def summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if force:
        prefix = summary_path + 'Force_' + model + '_'
    else:
        prefix = summary_path + model + '_'

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
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'MG94_clock_summary.txt',
                     unfinished_list_file = summary_path + 'MG94_clock_unfinished.txt',
                     clock = True, force = False)

    # MG94 nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'MG94_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'MG94_nonclock_unfinished.txt',
                     clock = False, force = False)

    # MG94 clock model (force tau = 0)
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_MG94_clock_summary.txt',
                     unfinished_list_file = summary_path + 'Force_MG94_clock_unfinished.txt',
                     clock = True, force = True)
    # MG94 nonclock model (force tau = 0)
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_MG94_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'Force_MG94_nonclock_unfinished.txt',
                     clock = False, force = True)
                
####################################################################################################################################################
    model = 'HKY'
    # HKY clock model
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'HKY_clock_summary.txt',
                     unfinished_list_file = summary_path + 'HKY_clock_unfinished.txt',
                     clock = True, force = False)

    # HKY nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'HKY_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'HKY_nonclock_unfinished.txt',
                     clock = False, force = False)

    # HKY clock model (force tau = 0)
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_HKY_clock_summary.txt',
                     unfinished_list_file = summary_path + 'Force_HKY_clock_unfinished.txt',
                     clock = True, force = True)
    # HKY nonclock model (force tau = 0)
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_HKY_nonclock_summary.txt',
                     unfinished_list_file = summary_path + 'Force_HKY_nonclock_unfinished.txt',
                     clock = False, force = True)
                
