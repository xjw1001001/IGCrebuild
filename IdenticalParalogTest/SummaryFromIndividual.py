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


    summary_path = '/Users/xji3/IdenticalParalogTestFromCluster10072015/Summary/'
    
    model = 'MG94'
    clock = False
    force = False
    summary_file = summary_path + model + '_nonclock_summary.txt'
    unfinished_list_file = summary_path + model + '_nonclock_unfinished.txt'
    
    summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force)

        

        
