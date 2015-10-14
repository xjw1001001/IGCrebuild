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


def summary_from_ind(paralog, BootNumList, summary_path, model, clock, force):
    summary_file = gen_summary_file_name(paralog, summary_path, model, clock, force)
    summary_mat = []
    
    for boot_num in BootNumList:
        ind_summary_file = summary_path + '_'.join(paralog) + '_MG94_nonclock_Sim_' + str(boot_num) + '_summary.txt'
        if os.path.isfile(ind_summary_file):
            res = np.loadtxt(open(ind_summary_file, 'r'))
            if len(np.atleast_1d(res)) > 1:
                summary_mat.append(res.tolist())
                label = open(ind_summary_file, 'r').readlines()[-1][2:-1]

    t = np.matrix(summary_mat)
    header = ' '.join(['Simulation' + str(boot_num) for boot_num in BootNumList])
    footer = label
    np.savetxt(open(summary_file, 'w+'), t.T, delimiter = ' ', header = header, footer = footer)


if __name__ == '__main__':

    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
            
    #pairs = [pairs[0]]
    for paralog in pairs:
        summary_path = '/Users/xji3/SimulationFromCluster10122015/SimulationSummary/' + '_'.join(paralog) + '/'
        model = 'MG94'
        clock = False
        force = False
        BootNumList = range(1, 101)
        
        summary_from_ind(paralog, BootNumList, summary_path, model, clock, force)
    
        

        
