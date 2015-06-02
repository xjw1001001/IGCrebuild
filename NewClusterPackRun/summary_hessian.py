import numpy as np
import os

def summary_from_ind(pairs, summary_path, model, summary_file, unfinished_list_file, clock, force, Dir, gBGC):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    edge_list = ['N0_N1', 'N0_kluyveri', 'N1_N2', 'N1_castellii', 'N2_N3', 'N2_bayanus', 'N3_N4', 'N3_kudriavzevii', 'N4_N5', 'N4_mikatae', 'N5_cerevisiae', 'N5_paradoxus']
    process_list = ['%AG', '%A', '%C', 'kappa']
    if model == 'MG94':
        process_list.append('omega')
        
    if force:
        prefix = summary_path + 'Force_' + model + '_'
    else:
        if Dir:
            if gBGC:
                prefix = summary_path + 'gBGC_Dir_' + model + '_'
                process_list.extend(['tau12', 'tau21', 'gamma'])
            else:
                prefix = summary_path + 'Dir_' + model + '_'
                process_list.extend(['tau12', 'tau21'])
        else:
            if gBGC:
                prefix = summary_path + 'gBGC_' + model + '_'
                process_list.extend(['tau', 'gamma'])
            else:
                prefix = summary_path + model + '_'
                process_list.append('tau')

    if clock:
        suffix = '_clock_hessian.txt'
    else:
        suffix = '_nonclock_hessian.txt'

    para_list = np.concatenate((process_list, edge_list))
    label = ' '.join(para_list)
    for i in range(len(para_list)):
        for j in range(i + 1, len(para_list)):
            label += ' ' + para_list[i] + ',' + para_list[j]
    for pair in pairs:
        var_summary = []
        summary = prefix + '_'.join(pair) + suffix
        if os.path.isfile(summary):
            res = np.loadtxt(open(summary, 'r'))
            assert(res.shape[0] == len(para_list))
            var_mat = np.linalg.inv( - res)
            var_summary.extend(var_mat.diagonal())
            for i in range(len(para_list)):
                for j in range(i + 1, len(para_list)):
                    var_summary.append(var_mat[i, j])
                
            if len(np.atleast_1d(res)) > 1:
                summary_mat.append(var_summary)
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
    summary_path = '/Users/xji3/FromCluster06012015/NewPackageNewRun/'
#    summary_path = '/Users/xji3/FromCluster05082015/NewPackageNewRun04102015/'
    
    pairs = []
    all_pairs = '../All_Pairs.txt'
    jeff_pairs = './Jeff_pairs_list.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])

    #pairs = [pairs[0]]
    model = 'HKY'
    summary_file = summary_path + 'HKY_nonclock_hessian_summary.txt'
    unfinished_list_file = summary_path + 'HKY_nonclock_hessian_unfinished.txt'
    clock = False
    force = False
    Dir = False
    gBGC = False

    model = 'MG94'

    # MG94 nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'MG94_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'MG94_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = False, gBGC = False)

    # MG94 Dir nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_MG94_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'Dir_MG94_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = True, gBGC = False)

    # MG94 gBGC nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_MG94_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'gBGC_MG94_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = False, gBGC = True)

    # MG94 Dir gBGC nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_MG94_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_MG94_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = True, gBGC = True)

##    # MG94 nonclock model (force tau = 0)
##    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_MG94_nonclock_hessian.txt',
##                     unfinished_list_file = summary_path + 'Force_MG94_nonclock_hessian_unfinished.txt',
##                     clock = False, force = True, Dir = False, gBGC = False)
                
####################################################################################################################################################
    model = 'HKY'

    # HKY nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'HKY_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'HKY_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = False, gBGC = False)
    # HKY Dir nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_HKY_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'Dir_HKY_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = True, gBGC = False)

    # HKY gBGC nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'gBGC_HKY_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'gBGC_HKY_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = False, gBGC = True)

    # HKY Dir gBGC nonclock model 
    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Dir_gBGC_HKY_nonclock_hessian.txt',
                     unfinished_list_file = summary_path + 'Dir_gBGC_HKY_nonclock_hessian_unfinished.txt',
                     clock = False, force = False, Dir = True, gBGC = True)

##    # HKY nonclock model (force tau = 0)
##    summary_from_ind(pairs, summary_path, model, summary_file = summary_path + 'Force_HKY_nonclock_hessian.txt',
##                     unfinished_list_file = summary_path + 'Force_HKY_nonclock_hessian_unfinished.txt',
##                     clock = False, force = True, Dir = False, gBGC = False)

    sh_line = 'sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/'
    with open('./unfinished_hessian.sh', 'w+') as f:
        f.write('#!/bin/bash' + '\n')
        for model in ['HKY', 'MG94']:
            for case in ['', 'Dir', 'gBGC', 'Dir_gBGC']:
                if case == '':
                    unfinished_list_file = summary_path + case + model + '_nonclock_hessian_unfinished.txt'
                else:
                    unfinished_list_file = summary_path + case + '_' + model + '_nonclock_hessian_unfinished.txt'
                if 'Dir' in case:
                    Dir = '--dir'
                else:
                    Dir = '--no-dir'

                if 'gBGC' in case:
                    gBGC = '--gBGC'
                else:
                    gBGC = '--no-gBGC'

                clock = '--no-clock'
                if 'Force' in case:
                    force = '--force'
                else:
                    force = '--no-force'

                unfinished_pairs = []
                
                with open(unfinished_list_file, 'r') as s:
                    for line in s.readlines():
                        unfinished_pairs.append(line.replace('\n','').split('_'))

                for pair in unfinished_pairs:
                    f.write(sh_line + '_'.join(pair) + '_'.join(['', case, model, 'Hessian_unfinished.sh']) + '\n')
                    with open('./NewRun/' + '_'.join(pair) + '_'.join(['', case, model, 'Hessian_unfinished.sh']),'w+') as g:
                        g.write('#!/bin/bash' + '\n')
                        g.write('python GenerateInfomationMatrix.py --paralog1 ' + pair[0] + ' --paralog2 ' + pair[1] + ' '.join(['', '--model', model, clock, force, Dir, gBGC]) + '\n')
                    
                                
