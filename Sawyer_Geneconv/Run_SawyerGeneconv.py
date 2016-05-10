import os
import subprocess

def initialize(paralog, out_path = './output/', alignment_path = '../MafftAlignment/'):
    if not os.path.isdir(out_path + '_'.join(paralog)):
        os.mkdir(out_path + '_'.join(paralog))

    input_alignment = out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta' 
    if not os.path.isfile(input_alignment):
        subprocess.check_output(['cp', input_alignment.replace(out_path, alignment_path), input_alignment]) 

def run_SawyerGeneconv(paralog, out_path = './output/', software_dir = '/Users/xji3/Genconv/Sawyer_Geneconv/geneconv', protein = False):
    if protein:
        SawyerGeneconvCmd = [software_dir,out_path + '_'.join(paralog) + '/' +  '_'.join(paralog) + '_input.fasta', '-nolog', '/lp', '/w822', '/p']
    else:
        SawyerGeneconvCmd = [software_dir,out_path + '_'.join(paralog) + '/' +  '_'.join(paralog) + '_input.fasta', '-nolog', '/lp', '/w822']
    subprocess.check_output(SawyerGeneconvCmd)

def read_SawyerGeneconv(paralog, file_path = './output/'):
    file_name = file_path + '_'.join(paralog) + '/' +  '_'.join(paralog) + '_input.frags'
    with open(file_name, 'r') as f:
        content = f.readlines()

    pline = 0
    for i in range(2, len(content)):
        #print content[i][:7], content[i][:7] == '# SCORE', content[i - 1][:8] == '#  frags', content[i - 2][:8] == '# Inner '
        if content[i][:7] == '# SCORE' and content[i - 1][:8] == '#  frags' and content[i - 2][:8] == '# Inner ':
            pline = i
    return float(content[pline].split()[3])

if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))

    
    for paralog in pairs:
        initialize(paralog)
        run_SawyerGeneconv(paralog)
        print '_'.join(paralog), read_SawyerGeneconv(paralog)
        run_SawyerGeneconv(paralog, protein = True)
        print '_'.join(paralog), read_SawyerGeneconv(paralog)
        

