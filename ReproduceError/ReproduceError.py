from Rewrite_CodonGeneconv import ReCodonGeneconv
import numpy as np
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    #num_boot = args.bootnum
    Force = None
    clock = False
    alignment_file = '../MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    save_name = None
    newicktree = './YeastTree.newick'
    omega_guess = 0.1

    
    test_hky = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = False, save_name = save_name)
    test_hky.get_mle(False, True, 0, 'BFGS')

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = False)
    x = np.concatenate((test_hky.x_process[:-1], np.log([omega_guess]), test_hky.x_process[-1:], test_hky.x_rates))
    test.update_by_x(x)
    test.get_mle(True, True, 0, 'BFGS')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    #parser.add_argument('--bootnum', type = int, help = 'Number of Bootstrap Replicate')
    
    main(parser.parse_args())
##    
##    pairs = []
##    all_pairs = './Filtered_pairs.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##    
##    model = 'MG94_'
##    sh_line = 'sbatch -o Boot-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
##
##    pairs = [pairs[0]]
##    
##    with open('./Run_Boot.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for paralog in pairs:
##            for num_boot in range(100):
##                f.write(sh_line + model + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '.sh\n')
##                with open('./ShFiles/' + model + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    printscreen = ' > ' + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '_PrintScreen.txt\n'
##                    g.write('python RunBootstrap.py ' + ' --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --bootnum ' + str(num_boot + 1) + printscreen)
##
##

    
    
