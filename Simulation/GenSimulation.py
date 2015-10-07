from GeneconvSimulation import *
import os

if __name__ == '__main__':

    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    newicktree = './YeastTree.newick'

    pairs = [pairs[0]]
    for paralog in pairs:
        if not os.path.isdir('./Simulation/' + '_'.join(paralog)):
            os.mkdir('./Simulation/' + '_'.join(paralog))
            
        if not os.path.isdir('./SimulationSave/' + '_'.join(paralog)):
            os.mkdir('./SimulationSave/' + '_'.join(paralog))

        if not os.path.isdir('./SimulationSummary/' + '_'.join(paralog)):
            os.mkdir('./SimulationSummary/' + '_'.join(paralog))

        x = np.exp(np.loadtxt(open('./save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt', 'r')))
        for sim_num in range(10, 101):
            test = SimGeneconv(newicktree, paralog, x, Model = 'MG94', nnsites = 112, Dir = False, gBGC = False)
            test.sim()
            
            test.output_seq(path = './Simulation/' + '_'.join(paralog) + '/', sim_num = sim_num)
            
            with open('./Simulation/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_'.join([test.Model, 'non_clock', 'sim', str(sim_num), 'TrueEvent.txt']), 'w+') as f:
                f.write('\t'.join(['edge', '# IGC', '# All', '\n']))
                for edge in test.edge_list:
                    f.write('\t'.join(['_'.join(edge), str(test.node_to_sim[edge[1]][1]), str(test.node_to_sim[edge[1]][2])]) + '\n')
                actual_ratio = sum([test.node_to_sim[i][1] + 0.0 for i in test.node_to_sim if i != 'kluyveri']) / sum([test.node_to_sim[i][2] + 0.0 for i in test.node_to_sim if i != 'kluyveri'])
                f.write(str(actual_ratio) + '\n')

