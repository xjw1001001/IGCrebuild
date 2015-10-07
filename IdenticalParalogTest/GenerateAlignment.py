import os
import subprocess
from Bio import SeqIO

if __name__ == '__main__':
    in_path = '/Users/xji3/Genconv/MafftAlignment/'
    out_path = '/Users/xji3/Genconv/IdenticalParalogAlignment/'
    pairs = []
    with open('./All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    if ['YLR028C', 'YMR120C'] in pairs:
        pairs.remove(['YLR028C', 'YMR120C'])

    #pairs = [pairs[0]]
    for paralog in pairs:
        if not os.path.isdir(out_path + '_'.join(paralog)):
            os.mkdir(out_path + '_'.join(paralog))

        in_alignment = in_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
        out_alignment = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_IdenticalParalog_input.fasta'
        paml_out_alignment = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_IdenticalParalog_paml_input.fasta'

        
        seq_dict = SeqIO.to_dict(SeqIO.parse( in_alignment, "fasta" ))
        name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}

        for name in name_to_seq.keys():
            if name[-len(paralog[0]):] == paralog[1]:
                other_paralog_name = name[:-len(paralog[1])] + paralog[0]
                name_to_seq[name] = name_to_seq[other_paralog_name]
                print name + ' sequence overwrited by ' + other_paralog_name

        with open(out_alignment, 'w+') as f:
            for name in name_to_seq.keys():
                f.write('>' + name + '\n')
                f.write(name_to_seq[name] + '\n')
                
        with open(paml_out_alignment, 'w+') as f:
            for name in name_to_seq.keys():
                if name[-len(paralog[0]):] == paralog[0]:
                    f.write('>' + name + '\n')
                    f.write(name_to_seq[name] + '\n')
                
