'''
New ortholog mapping script
'''
from YeastModel import *
import subprocess


def read_ortholog(input_file_dir = '/Users/xji3/blast/YeastAlignment/all-output.txt'):
    one_to_one_dict = {}
    many_to_many_dict = {}
    with open(input_file_dir,'r') as f:
        read_now = False
        for line in f:
            if line.split()[0] =='>':
                read_now = True
                continue
            if read_now:
                info = line.split()
                if info[1]==info[2]:
                    if info[3][0:4]=='Scer':
                        one_to_one_dict[info[3][5:]]={info[i][0:4]:info[i][5:] for i in range(4,int(info[1]))}
                else:
                    if info[3][0:4]=='Scer':
                        key = [info[3][5:]]
                        j = 4
                        while info[j][0:4]=='Scer':
                            key.append(info[j][5:])
                            j += 1
                        many_to_many_dict[tuple(set(key))] = {}
                        for i in range( 3+len(key),int(info[1])):
                            if many_to_many_dict[tuple(set(key))].has_key(info[i][0:4]):
                                many_to_many_dict[tuple(set(key))][info[i][0:4]].append(info[i][5:])
                            else:
                                many_to_many_dict[tuple(set(key))][info[i][0:4]]=[info[i][5:]]
                            
                read_now = False

    return one_to_one_dict,many_to_many_dict

def processAlignment(input_file):
    align = AlignIO.read(input_file,'fasta')
    i=0
    while(i<align.get_alignment_length()):
        if not align[:,i].find('-') == -1:
            if i==0:
                align = align[:,1:]
            elif i==(align.get_alignment_length()-1):
                align = align[:,:-1]
            else:
                align = align[:,:i]+align[:,(i+1):]
        else:
            i=i+1
    assert(align.get_alignment_length()%3==0)
    return align

            

if __name__=='__main__':
    pair_list_dir = './data/Pairs/pairs_from_paper.txt'
    species_list = ['bayanus','cerevisiae','kudriavzevii','mikatae','paradoxus','castellii','kluyveri']
    gene_pairs = readInitialPairs()
    pair_list = gene_pairs.keys()

###
#Read ortholog info from http://www.broadinstitute.org/regev/orthogroups/
###
    one_to_one_dict,many_to_many_dict = read_ortholog()

###
#Read orthlog mapping for 5 species from Hittinger paper 'The awesome power of Yeast evolutionary genetics' supplementary material
###
    wb = open_workbook('/Users/xji3/blast/YeastAlignment/HittingerSequence/Supp_Table1_OrthologSets.xls')
    sheet = wb.sheets()[0]
    Gene2map = {}
    col_list = [str(s.value) for s in sheet.row(0)]
    for row in range(1,sheet.nrows):
        gene_name = str(sheet.cell(row,4).value)
        info = {col_list[col]:sheet.cell(row,col).value for col in range(5,sheet.ncols)}
        Gene2map[gene_name]=info

###
#Intersection between two mappings
###
    filtered_one_to_one_dict = {key:one_to_one_dict[key] for key in one_to_one_dict.keys() if one_to_one_dict[key].has_key('Scas') and one_to_one_dict[key].has_key('Sklu') }
    filtered_many_to_many_dict = {key:many_to_many_dict[key] for key in many_to_many_dict.keys() if many_to_many_dict[key].has_key('Scas') and many_to_many_dict[key].has_key('Sklu') }



    filtered_pairs = []
    for pair in pair_list:
        paralog1 = pair[0]
        paralog2 = pair[1]
        if paralog1 in Gene2map.keys() and paralog2 in Gene2map.keys():
            filtered_pairs.append(pair)

    filtered_pair_list = [key for key in filtered_pairs if key in many_to_many_dict.keys() or (key[1],key[0]) in many_to_many_dict.keys()]
    for pair in filtered_pair_list:
        if not pair in many_to_many_dict:
            filtered_pair_list.remove(pair)
            filtered_pair_list.append((pair[1],pair[0]))

    Filtered_one_outgroup = [pair for pair in filtered_pair_list if pair in filtered_many_to_many_dict and len(filtered_many_to_many_dict[pair]['Sklu']) == 1 and  len(filtered_many_to_many_dict[pair]['Scas']) == 2]
    Filtered_two_outgroup = [pair for pair in filtered_pair_list if pair in filtered_many_to_many_dict and len(filtered_many_to_many_dict[pair]['Sklu']) == 1 and  len(filtered_many_to_many_dict[pair]['Scas']) == 1]

###
#Extract DNA sequences of each species for pairs
###
    castellii_data = {}
    castellii_name2id = {}
    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/NT.fsa'
    handle = open(input_file,'rU')
    for record in SeqIO.parse(handle,'fasta'):
        if 'NCAS' in record.id:
            castellii_data[record.id]=str(record.seq)
            for item in record.description.split():
                if item[0] == 'Y':
                    castellii_name2id[item]=record.id
                    break
                
    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/Ncastellii_genome.tab'
    with open(input_file,'rU') as handle:
        for line in handle:
            if line.split()[-1][0] =='Y':
                if not line.split()[-1] in castellii_name2id:
                    castellii_name2id[line.split()[-1]]=line.split()[0]
                    print line.split()[-1]

    Pair_to_sequence = {pair:{pair[0]:{},pair[1]:{}} for pair in Filtered_one_outgroup if pair[0] in castellii_name2id or pair[1] in castellii_name2id}
    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/Scas.fasta'
    data_castellii = {}
    handle = open(input_file,'rU')
    for record in SeqIO.parse(handle,'fasta'):
        data_castellii[record.id]=str(record.seq)

    remove_pair_list = []

    for pair in Pair_to_sequence:
        if pair[0] in castellii_name2id:
            seq_known = castellii_data[castellii_name2id[pair[0]]]
            Pair_to_sequence[pair][pair[0]]['castellii']=seq_known
            [id1,id2]=filtered_many_to_many_dict[pair]['Scas']
            if data_castellii[id1] == seq_known:
                Pair_to_sequence[pair][pair[1]]['castellii']=data_castellii[id2]
            elif data_castellii[id2] == seq_known:
                Pair_to_sequence[pair][pair[1]]['castellii']=data_castellii[id1]
            else:
                print '=====================Warning====================='
                print pair,' cannot find theirs sequences'
                remove_pair_list.append(pair)
        else:
            seq_known = castellii_data[castellii_name2id[pair[1]]]
            Pair_to_sequence[pair][pair[1]]['castellii']=seq_known
            [id1,id2]=filtered_many_to_many_dict[pair]['Scas']
            if data_castellii[id1] == seq_known:
                Pair_to_sequence[pair][pair[0]]['castellii']=data_castellii[id2]
            elif data_castellii[id2] == seq_known:
                Pair_to_sequence[pair][pair[0]]['castellii']=data_castellii[id1]
            else:
                print '=====================Warning====================='
                print pair,' cannot find theirs sequences'
                remove_pair_list.append(pair)            
    for pair in remove_pair_list:
        Pair_to_sequence.pop(pair)

    species_list.remove('castellii')
    species_list.remove('kluyveri')
    data_5species = {spe:{} for spe in species_list}
    for spe in species_list:
        with open('/Users/xji3/blast/YeastAlignment/HittingerSequence/S'+spe[0:3]+'.fsa','rU') as f:
            for record in SeqIO.parse(f,'fasta'):
                data_5species[spe][record.id] = str(record.seq)
    for pair in Pair_to_sequence:
        for paralog in pair:
            for spe in species_list:
                Pair_to_sequence[pair][paralog][spe] = data_5species[spe][str(Gene2map[paralog]['S'+spe[0:3]+'_Name'])]


    data_kluyveri = {}
    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/Sklu.fasta'
    with open(input_file,'rU') as handle:
        for record in SeqIO.parse(handle,'fasta'):
            data_kluyveri[record.id] = str(record.seq)

    for pair in Pair_to_sequence:
        spe = 'kluyveri'
        Pair_to_sequence[pair][spe] = data_kluyveri[filtered_many_to_many_dict[pair]['Sklu'][0]]

###
#Now filter by sequence length before making alignment
###
    species_list = ['bayanus','cerevisiae','kudriavzevii','mikatae','paradoxus','castellii']
    count = 0
    remove_pair_list = []
    for pair in Pair_to_sequence:
        for paralog in pair:
            length_list = [len(Pair_to_sequence[pair][paralog][spe]) for spe in species_list]
            length_list.append(len(Pair_to_sequence[pair]['kluyveri']))
            max_length = max(length_list)+0.0
            min_length = min(length_list)+0.0
            if min_length/max_length < 0.9:
                if not pair in remove_pair_list:
                    remove_pair_list.append(pair)
                count += 1
    print 'Out of ',len(Pair_to_sequence),'pairs, '+str(len(remove_pair_list))+' pairs removed'

    Pairs_for_alignment = deepcopy(Pair_to_sequence.keys())
    for pair in remove_pair_list:
        Pairs_for_alignment.remove(pair)

###
#Now output fasta file for alignment input
###
    for pair in Pairs_for_alignment:
        paralog1 = pair[0]
        paralog2 = pair[1]
        if not os.path.exists('./PairsAlignemt/'+paralog1+'_'+paralog2+'/'):
            os.mkdir('./PairsAlignemt/'+paralog1+'_'+paralog2+'/')
        with open('./PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'.fa','w+') as f:
            for spe in species_list:
                f.write('>'+spe+paralog1+'\n')
                f.write(Pair_to_sequence[pair][paralog1][spe]+'\n')
                f.write('>'+spe+paralog2+'\n')
                f.write(Pair_to_sequence[pair][paralog2][spe]+'\n')
            f.write('>'+'kluyveri'+paralog1+'\n')
            f.write(Pair_to_sequence[pair]['kluyveri']+'\n')
    
###
#Now use fsa to do multiple sequence alignment
###
    for pair in Pairs_for_alignment:
        paralog1 = pair[0]
        paralog2 = pair[1]
        #cmd = ['/Users/xji3/fsa-1.15.9/src/main/fsa','/Users/xji3/Genconv/PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'.fa','>','/Users/xji3/Genconv/PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_fsa_alignment.fa']
        #cmd = ['fsa','/Users/xji3/Genconv/PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'.fa','>','/Users/xji3/Genconv/PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_fsa_alignment.fa']
        #subprocess.check_output(cmd)


        cmd = ['perl','./translatorx_vLocal.pl','-i','./PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'.fa','-o','./PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_translatorX_alignment','-p','M']
        #os.system(' '.join(cmd))
        my_env = os.environ.copy()
        my_env['PATH'] = "/usr/sbin:/sbin:" + '/Users/xji3/AlignmentPro:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/xji3/fas-1.15.9/MUMmer3.23:/opt/X11/bin:/opt/ImageMagick/bin:/usr/local/ncbi/blast/bin:/usr/texbin:/usr/local/ncbi/blast/bin:/Users/xji3/AlignmentPro'
        subprocess.check_output(cmd,env = my_env)

        aligned_file = './PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_translatorX_alignment.nt_ali.fasta'
        align = processAlignment(aligned_file)
        lookup = dict((rec.id,str(rec.seq)) for rec in align)
        with open('./PairsAlignemt/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_input.fasta','w+') as f:
            for rec in align:
                f.write('>'+str(rec.id)+'\n')
                f.write(str(rec.seq)[:-3]+'\n')
###
#Now generate python code for cluster
###
    Sample_py_rooted_free = './PairsAlignemt/Rooted_Codon_Free_Tau_YBL087C_YER117W.py'
    Sample_py_rooted_force = './PairsAlignemt/Rooted_Codon_Force_Tau_YBL087C_YER117W.py'
    Sample_py_unrooted_free = './PairsAlignemt/UnRooted_Codon_Free_Tau_YBL087C_YER117W.py'
    Sample_py_unrooted_force = './PairsAlignemt/UnRooted_Codon_Force_Tau_YBL087C_YER117W.py'

    replace_paralog1 = 'YBL087C'
    replace_paralog2 = 'YER117W'

    for pair in Pairs_for_alignment:
        paralog1 = pair[0]
        paralog2 = pair[1]
        new_py_rooted_free = Sample_py_rooted_free.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt','./PairsAlignemt/'+paralog1+'_'+paralog2)
        new_py_rooted_force = Sample_py_rooted_force.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt','./PairsAlignemt/'+paralog1+'_'+paralog2)
        new_py_unrooted_free = Sample_py_unrooted_free.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt','./PairsAlignemt/'+paralog1+'_'+paralog2)
        new_py_unrooted_force = Sample_py_unrooted_force.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt','./PairsAlignemt/'+paralog1+'_'+paralog2)

        with open(Sample_py_rooted_free,'rU') as f:
            with open(new_py_rooted_free,'w+') as g:
                for line in f:
                    g.write(line.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2))

        with open(Sample_py_rooted_force,'rU') as f:
            with open(new_py_rooted_force,'w+') as g:
                for line in f:
                    g.write(line.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2))
        
        with open(Sample_py_unrooted_free,'rU') as f:
            with open(new_py_unrooted_free,'w+') as g:
                for line in f:
                    g.write(line.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2))

        with open(Sample_py_unrooted_force,'rU') as f:
            with open(new_py_unrooted_force,'w+') as g:
                for line in f:
                    g.write(line.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2))

###
#Now generate sh file for cluster
###
    with open('./PairsAlignemt/run_all_yeast_pairs.sh','w+') as f:
        f.write('#!/bin/bash'+'\n')
        for pair in Pairs_for_alignment:
            paralog1 = pair[0]
            paralog2 = pair[1]
            py_list = [Sample_py_rooted_free,Sample_py_rooted_force,Sample_py_unrooted_free,Sample_py_unrooted_force]
            for py_file in py_list:
                target_py_file = py_file.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt/','')
                sh_file = py_file.replace(replace_paralog1,paralog1).replace(replace_paralog2,paralog2).replace('./PairsAlignemt/','./PairsAlignemt/'+paralog1+'_'+paralog2+'/Run_').replace('.py','_script.sh')
                f.write(' '.join(['sbatch -o CodonGeneconv-%j.out -p bigmem -w node92 --mail-type=FAIL --mail-user=xji3@ncsu.edu',sh_file.replace('./PairsAlignemt/','./')])+'\n')
                with open(sh_file,'w+') as g:
                    g.write('#!/bin/bash'+'\n')
                    g.write('cd '+paralog1+'_'+paralog2+'\n')
                    g.write(' '.join(['python',target_py_file,'>',target_py_file.replace('.py','_PrintScreen.txt')]))
