from GeneconvSimulation import *
import os
from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import numpy as np
import random
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]

    pairs = []
    pairs_to_nsites = {}
    all_pairs = './Filtered_pairs_length.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            info = line.replace('\n','').split('_')
            pairs.append(info[:2])
            pairs_to_nsites['_'.join(info[:2])] = int(info[-1])

    newicktree = './YeastTree.newick'
    np.random.seed(822)
    random.seed(822)
    
    if not os.path.isdir('./Simulation/' + '_'.join(paralog)):
        os.mkdir('./Simulation/' + '_'.join(paralog))
        
    if not os.path.isdir('./SimulationSave/' + '_'.join(paralog)):
        os.mkdir('./SimulationSave/' + '_'.join(paralog))

    if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog)):
        os.mkdir('./SimulationSummary/' + '_'.join(paralog))

    x = np.exp(np.loadtxt(open('./save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt', 'r')))
    for sim_num in range(1, 101):
        test = SimGeneconv(newicktree, paralog, x, Model = 'MG94', nnsites = pairs_to_nsites['_'.join(paralog)], Dir = False, gBGC = False)
        test.sim()
        
        test.output_seq(path = './Simulation/' + '_'.join(paralog) + '/', sim_num = sim_num)
        
        with open('./Simulation/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_'.join([test.Model, 'non_clock', 'sim', str(sim_num), 'TrueEvent.txt']), 'w+') as f:
            f.write('\t'.join(['edge', '# IGC', '# All', '\n']))
            for edge in test.edge_list:
                f.write('\t'.join(['_'.join(edge), str(test.node_to_sim[edge[1]][1]), str(test.node_to_sim[edge[1]][2])]) + '\n')
            actual_ratio = sum([test.node_to_sim[i][1] + 0.0 for i in test.node_to_sim if i != 'kluyveri']) / sum([test.node_to_sim[i][2] + 0.0 for i in test.node_to_sim if i != 'kluyveri'])
            f.write(str(actual_ratio) + '\n')

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    
    main(parser.parse_args())
## 
##    pairs = []
##    pairs_to_nsites = {}
##    all_pairs = './Filtered_pairs_length.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            info = line.replace('\n','').split('_')
##            pairs.append(info[:2])
##            pairs_to_nsites['_'.join(info[:2])] = int(info[-1])
##
##    newicktree = './YeastTree.newick'
##    np.random.seed(822)
##    random.seed(822)
##
##    for paralog in pairs:
##        if not os.path.isdir('./Simulation/' + '_'.join(paralog)):
##            os.mkdir('./Simulation/' + '_'.join(paralog))
##            
##        if not os.path.isdir('./SimulationSave/' + '_'.join(paralog)):
##            os.mkdir('./SimulationSave/' + '_'.join(paralog))
##
##        if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog)):
##            os.mkdir('./SimulationSummary/' + '_'.join(paralog))
##
##        x = np.exp(np.loadtxt(open('./save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt', 'r')))
##        for sim_num in range(1, 101):
##            test = SimGeneconv(newicktree, paralog, x, Model = 'MG94', nnsites = pairs_to_nsites['_'.join(paralog)], Dir = False, gBGC = False)
##            test.sim()
##            
##            test.output_seq(path = './Simulation/' + '_'.join(paralog) + '/', sim_num = sim_num)
##            
##            with open('./Simulation/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_'.join([test.Model, 'non_clock', 'sim', str(sim_num), 'TrueEvent.txt']), 'w+') as f:
##                f.write('\t'.join(['edge', '# IGC', '# All', '\n']))
##                for edge in test.edge_list:
##                    f.write('\t'.join(['_'.join(edge), str(test.node_to_sim[edge[1]][1]), str(test.node_to_sim[edge[1]][2])]) + '\n')
##                actual_ratio = sum([test.node_to_sim[i][1] + 0.0 for i in test.node_to_sim if i != 'kluyveri']) / sum([test.node_to_sim[i][2] + 0.0 for i in test.node_to_sim if i != 'kluyveri'])
##                f.write(str(actual_ratio) + '\n')
##
