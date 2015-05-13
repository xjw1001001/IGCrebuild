from Rewrite_CodonGeneconv import ReCodonGeneconv
from GeneconvSimulation import SimGeneconv
import numpy as np

if __name__ == '__main__':

    paralog1 = 'YDR502C'
    paralog2 = 'YLR180W'

    paralog = [paralog1, paralog2]
    newicktree = '../PairsAlignemt/YeastTree.newick'
    log_omega = np.log(8.562946237878546474e-02)
    x = np.exp([-0.76043136, -0.56202514, -0.86577248,  1.34559939, -3.18120517,
       -0.53628023, -1.58482627, -1.45313697, -1.48522716, -0.51722434,
       -2.40816466, -2.06956473, -2.4618848 , -1.75893577, -2.61506075,
       -1.66772197, -2.38438532, -2.40841467])

    #YBR117C_YPR074C
    x_HKY = np.exp([-0.75572536, -0.57064656, -0.78114092,  1.43680162, -0.12676788,
       -1.79601116, -3.14039449, -2.25513103, -1.79522678, -3.6940989 ,
       -3.04007122, -3.73703586, -3.06174264, -4.02058201, -2.91143155,
       -3.28320901, -3.56320802])

    x_HKY[4] = 0.0

    x_gBGC_dir = np.exp([-0.76332984, -0.6035293 , -0.77587976,  1.36553822, -3.2183019 ,
       -1.85738889, -0.22358248, -3.41664095, -1.54763648, -1.54711367,
       -1.43821404, -0.42501841, -2.32832516, -1.94156997, -2.37697147,
       -1.7136708 , -2.50543922, -1.64307197, -2.38638795, -2.33371383])

    #Dir = True
    Dir = False
    #gBGC = True
    gBGC = False
    #simulation = SimGeneconv( newicktree, paralog, x_gBGC_dir, Model = 'MG94', nnsites = 381, Dir = Dir, gBGC = gBGC)
    simulation = SimGeneconv( newicktree, paralog, x_HKY, Model = 'HKY', nnsites = 3 * 381, Dir = Dir, gBGC = gBGC)
    #simulation = SimGeneconv( newicktree, paralog, x, Model = 'MG94', nnsites = 381, Dir = Dir, gBGC = gBGC)
    simulation.sim()
    simulation.output_seq()

    suffix = '_' + simulation.Model
    if Dir:
        suffix += '_dir'
    if gBGC:
        suffix += '_gBGC'
    alignment_file = './simulation/' + '_'.join(paralog) + '/' + '_'.join(paralog) + suffix + '.fasta'
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = False)
    test.get_mle(False, True, 1, 'BFGS')

    print simulation.x
    print np.exp(test.x)

    for edge in test.edge_list:
        print simulation.edge_to_blen[edge], test.edge_to_blen[edge], edge

    print simulation.pi, simulation.kappa, simulation.tau
    print test.pi, test.kappa, test.tau
    for node in simulation.node_to_sequence:
        print simulation.node_to_sequence[node][1].count('A') / 1143.0, simulation.node_to_sequence[node][1].count('C')/1143.0, simulation.node_to_sequence[node][1].count('G') / 1143.0, simulation.node_to_sequence[node][1].count('T')/1143.0

        self = simulation

