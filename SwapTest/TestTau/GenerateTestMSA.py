import os, subprocess
from Bio import Seq, SeqIO, AlignIO

if __name__ == '__main__':
    path = '/Users/xji3/Genconv/NewClusterPackRun/TestTau/'
    path_from = '/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/'
    pairs = []
    with open('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/pairs_list.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    pairs.remove(['YLR028C', 'YMR120C'])
#    pairs = [pairs[0]]
    for pair in pairs:

        switch_seq = {'cerevisiae' + pair[0]: 'cerevisiae' + pair[0],
                      'cerevisiae' + pair[1]: 'castellii' + pair[1],
                      'castellii' + pair[0] : 'castellii' + pair[0],
                      'castellii' + pair[1] : 'cerevisiae' + pair[1],
                      'kluyveri' + pair[0] : 'kluyveri' + pair[0] 
                      }

        align = AlignIO.read(path_from + '_'.join(pair) + '/' + '_'.join(pair) + '_input.fasta', 'fasta')
        name_to_seq = {str(rec.id) : str(rec.seq) for rec in align}
        mkdir_cmd = ['mkdir', path + '_'.join(pair)]
        os.system(' '.join(mkdir_cmd))
        with open(path + '_'.join(pair) + '/' + '_'.join(pair) + '_switched.fasta', 'w+') as f:
            for seq_id in switch_seq.keys():
                f.write('>' + seq_id + '\n')
                f.write(name_to_seq[switch_seq[seq_id]] + '\n')
            
    
