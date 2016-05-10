import numpy as np
import os

def gen_summary_file_name(paralog, summary_path, model, clock, force):
    prefix_summary = summary_path + model
    if force:
        prefix_summary = prefix_summary + '_Force'

    if clock:
        suffix_summary = '_clock_Simulation_summary.txt'
    else:
        suffix_summary = '_nonclock_Simulation_summary.txt'

    summary_file = prefix_summary + '_' + '_'.join(paralog) + suffix_summary
    return summary_file



def summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, swap = False):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    
    if force:
        prefix = 'Force_' + model + '_'
    else:
        prefix = model + '_'

    if swap:
        prefix = 'switched_' + prefix

    if clock:
        suffix = '_clock_summary.txt'
    else:
        suffix = '_nonclock_summary.txt'

    label = ''
    for pair in pairs:
        summary = summary_path + '_'.join(pair) + '/' + prefix + '_'.join(pair) + suffix
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

    

if __name__ == '__main__':

    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
            
    summary_path = '/Users/xji3/SwapTestFromCluster10122015/TestTau/'
    model = 'MG94'
    clock = False
    force = False
    summary_file = summary_path + 'MG94_nonclock_summary.txt'
    unfinished_list_file = summary_path + 'MG94_nonclock_unfinished.txt'
    swap = False
    
    summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, swap)
    
    force = True
    summary_file = summary_path + 'Force_MG94_nonclock_summary.txt'
    unfinished_list_file = summary_path + 'Force_MG94_nonclock_unfinished.txt'
    swap = False
    
    summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, swap)
            
    force = False
    summary_file = summary_path + 'switched_MG94_nonclock_summary.txt'
    unfinished_list_file = summary_path + 'switched_MG94_nonclock_unfinished.txt'
    swap = True
    
    summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, swap)

    force = True
    summary_file = summary_path + 'switched_Force_MG94_nonclock_summary.txt'
    unfinished_list_file = summary_path + 'switched_Force_MG94_nonclock_unfinished.txt'
    swap = True
    
    summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, swap)    

        
