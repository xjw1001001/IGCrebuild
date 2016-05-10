import os
import subprocess
from Bio import Seq, SeqIO, AlignIO
from Bio.Phylo.PAML import codeml, baseml
import numpy as np

if __name__ == '__main__':
    path = '/Users/xji3/Genconv/IdenticalParalogAlignment/'
    pairs = []
    with open('./All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])
    pairs.append(['YLR284C','YOR180C'])  # this pair didn't appear this time
    #pairs.remove(['YML026C', 'YDR450W'])# remove it for now


    #pairs = [pairs[-1]]

    tree_pair = ['YML026C', 'YDR450W']
    with open('./YeastTree_paml.newick', 'r') as f:
        all_tree_lines = f.readlines()

    with open('./codeml_tail.ctl', 'r') as f:
        all_codeml_ctl_lines = f.readlines()

    with open('./baseml_tail.ctl', 'r') as f:
        all_baseml_ctl_lines = f.readlines()
    

    codeml = '/Users/xji3/Downloads/paml4.8/bin/codeml'
    baseml = '/Users/xji3/Downloads/paml4.8/bin/baseml'
    
    for pair in pairs:
        print 'Now run paml on pair ' + ' '.join(pair)
        seqfile = path + '_'.join(pair) + '/' + '_'.join(pair) + '_IdenticalParalog_paml_input.fasta'
        treefile = path + '_'.join(pair) + '/' + '_'.join(pair) + '_tree.newick'
        with open(treefile, 'w+') as f:
            for line in all_tree_lines:
                new_line = line.replace(tree_pair[0], pair[0])
                new_line = new_line.replace(tree_pair[1], pair[1])
                f.write(new_line)

        outfile_codeml = path + '_'.join(pair) + '/' + '_'.join(pair) + '_IdenticalParalog_codeml'
        codeml_ctlfile = path + '_'.join(pair) + '/' + '_'.join(pair) + '_IdenticalParalog_codeml_control.ctl'
        with open(codeml_ctlfile, 'w+') as f:
            f.writelines(['seqfile = ' + seqfile + '\n', 'treefile = ' + treefile + '\n', 'outfile = ' + outfile_codeml + '\n'])
            f.writelines(all_codeml_ctl_lines)

        codeml_cmd = [codeml, '_'.join(pair) + '_IdenticalParalog_codeml_control.ctl']
        os.chdir(path + '_'.join(pair) + '/')
        #os.system(' '.join(codeml_cmd))
        subprocess.check_output(codeml_cmd)

        outfile_baseml = path + '_'.join(pair) + '/' + '_'.join(pair) + '_IdenticalParalog_baseml'
        baseml_ctlfile = path + '_'.join(pair) + '/' + '_'.join(pair) + '_IdenticalParalog_baseml_control.ctl'
        with open(baseml_ctlfile, 'w+') as f:
            f.writelines(['seqfile = ' + seqfile + '\n', 'treefile = ' + treefile + '\n', 'outfile = ' + outfile_baseml + '\n'])
            f.writelines(all_baseml_ctl_lines)

        baseml_cmd = [baseml, '_'.join(pair) + '_IdenticalParalog_baseml_control.ctl']
        subprocess.check_output(baseml_cmd)
        

##    summary_mat = []
##    finished_list = []
##    label = ['MG94_codeml_tree_length', 'MG94_codeml_lnL', 'MG94_codeml_omega', 'MG94_codeml_kappa',
##             'HKY_baseml_tree_length', 'HKY_baseml_lnL', 'HKY_baseml_kappa']
##    footer = ' '.join(label)
##    
##    for pair in pairs:
##        codeml_result = codeml.read('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/' + '_'.join(pair) + '/' + '_'.join(pair) + '_codeml')
##        baseml_result = baseml.read('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/' + '_'.join(pair) + '/' + '_'.join(pair) + '_baseml')
##        summary_mat.append([codeml_result['NSsites'][0]['tree length'],
##                            codeml_result['NSsites'][0]['lnL'],
##                            codeml_result['NSsites'][0]['parameters']['omega'],
##                            codeml_result['NSsites'][0]['parameters']['kappa'],
##                            baseml_result['tree length'],
##                            baseml_result['lnL'],
##                            baseml_result['parameters']['kappa']])
##        finished_list.append(pair)
##
##    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
##    np.savetxt(open('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/paml_summary.txt', 'w+'), np.matrix(summary_mat).T, delimiter = ' ', footer = footer, header = header)
