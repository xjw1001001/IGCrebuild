from IGCexpansion.CodonGeneconFunc import *
import subprocess

def BootSample(paralog, num_sample, in_path, out_path, model = 'MG94', save_file = None, save_out_path = None):
    if not os.path.isdir(out_path + '_'.join(paralog)):
        os.mkdir(out_path + '_'.join(paralog))

    seqloc = in_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    seq_dict = SeqIO.to_dict(SeqIO.parse( seqloc, "fasta" ))
    
    name_to_seq = {name:str(seq_dict[name].seq) for name in seq_dict.keys()}
    if model == 'MG94':
        for name in name_to_seq.keys():
            assert(len(name_to_seq[name]) % 3 == 0)
            tmp_seq = [name_to_seq[name][3 * j : 3 * j + 3] for j in range(len(name_to_seq[name]) / 3 )]
            name_to_seq[name] = tmp_seq
            
    nsites = len(name_to_seq[name])
    for i in range(num_sample):
        boot_sample = np.random.choice(range(nsites), size = (nsites), replace = True)
        boot_name_to_seq = {name:[name_to_seq[name][site_num] for site_num in boot_sample] for name in name_to_seq.keys()}
        with open(out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot' + str(i + 1) + '.fasta', 'w+') as f:
            for name in boot_name_to_seq.keys():
                f.write('>' + name + '\n')
                f.write(''.join(boot_name_to_seq[name]) + '\n')

        # copy save files if provided
        if os.path.isfile(save_file):
            if not os.path.isdir(save_out_path + '_'.join(paralog)):
                os.mkdir(save_out_path + '_'.join(paralog))
            subprocess.check_output(['cp', save_file, save_out_path + '_'.join(paralog) + '/' + '_'.join(paralog) + '_Boot' + str(i + 1) + '_save.txt'])
    


if __name__ == '__main__':
    pairs = []
    all_pairs = './Filtered_pairs.txt'
    with open(all_pairs, 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))


    num_sample = 100
    in_path = './MafftAlignment/'
    out_path = './BootStrapSamples/'

    for paralog in pairs:
        save_file = './save/MG94_' + '_'.join(paralog) + '_nonclock_save.txt'
        BootSample(paralog, num_sample, in_path, out_path, save_file = save_file, save_out_path = './BootStrapSave/')

        if not os.path.isdir('./BootStrapSummary/' + '_'.join(paralog)):
            os.mkdir('./BootStrapSummary/' + '_'.join(paralog))
    
