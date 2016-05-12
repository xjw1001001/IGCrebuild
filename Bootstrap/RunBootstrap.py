from IGCexpansion.CodonGeneconv import ReCodonGeneconv
import argparse

def main(args):
    paralog = [args.paralog1, args.paralog2]
    num_boot = args.bootnum
    Force = None
    clock = False
    alignment_file = './BootStrapSamples/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot' + str(num_boot) + '.fasta'
    save_name = './BootStrapSave/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot' + str(num_boot) + '_save.txt'
    newicktree = './YeastTree.newick'
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = None, clock = False, save_name = save_name)
    test.get_mle(False, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './BootStrapSummary/' + '_'.join(paralog) + '/', file_name = './BootStrapSummary/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot' + str(num_boot) + '_summary.txt')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
    parser.add_argument('--bootnum', type = int, help = 'Number of Bootstrap Replicate')
    
    main(parser.parse_args())
    
##    pairs = []
##    all_pairs = './Filtered_pairs.txt'
##    with open(all_pairs, 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##    
##    model = 'MG94_'
##    sh_line = 'sbatch -o Boot-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
##
##    #pairs = [pairs[0]]
##
##    for paralog in pairs:
##    
##        with open('./Run_Boot_'+ '_'.join(paralog) +'.sh', 'w+') as f:
##            f.write('#!/bin/bash' + '\n')
##            for num_boot in range(100):
##                f.write(sh_line + model + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '.sh\n')
##                with open('./ShFiles/' + model + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '.sh', 'w+') as g:
##                    g.write('#!/bin/bash' + '\n')
##                    printscreen = ' > ' + '_'.join(paralog) + '_Boot' + str(num_boot + 1) + '_PrintScreen.txt\n'
##                    g.write('python RunBootstrap.py ' + ' --paralog1 ' + paralog[0] + ' --paralog2 ' + paralog[1] + ' --bootnum ' + str(num_boot + 1) + printscreen)



    
    
