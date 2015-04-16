from YeastOrthologMapping import *

if __name__=='__main__':
    pair_list_dir = './data/Pairs/pairs_from_paper.txt'
    species_list = ['bayanus','cerevisiae','kudriavzevii','mikatae','paradoxus','castellii','kluyveri']
    gene_pairs = readInitialPairs()
    pair_list = gene_pairs.keys()

###
#Read ortholog info from http://www.broadinstitute.org/regev/orthogroups/
###
#One to One mapping : http://www.broadinstitute.org/regev/orthogroups//orthologs.html
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
    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/Scas.fasta'
    handle = open(input_file,'rU')
    for record in SeqIO.parse(handle,'fasta'):
        castellii_data[record.id]=str(record.seq)

    input_file = '/Users/xji3/blast/YeastAlignment/Broadinstitute/Scer-Scas-orthologs.txt'
    with open(input_file,'rU') as handle:
        for line in handle:
            info = line.split()
            castellii_name2id[info[0]] = info[1]

    Pair_to_sequence = {pair:{pair[0]:{},pair[1]:{}} for pair in Filtered_one_outgroup if pair[0] in castellii_name2id and pair[1] in castellii_name2id}

    remove_pair_list = []
    for pair in Pair_to_sequence:
        Pair_to_sequence[pair][pair[0]]['castellii'] = castellii_data[castellii_name2id[pair[0]]]
        Pair_to_sequence[pair][pair[1]]['castellii'] = castellii_data[castellii_name2id[pair[1]]]
    

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
    print 'Out of ',len(Pair_to_sequence),'pairs, '+str(len(remove_pair_list))+' pairs removed based on sequence length filter'

    Pairs_for_alignment = deepcopy(Pair_to_sequence.keys())
    for pair in remove_pair_list:
        Pairs_for_alignment.remove(pair)

###
#Now output fasta file for alignment input
#    output nex file for RevBayes input
###
    for pair in Pairs_for_alignment:
        paralog1 = pair[0]
        paralog2 = pair[1]
        if not os.path.exists('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/'+paralog1+'_'+paralog2+'/'):
            os.mkdir('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/'+paralog1+'_'+paralog2+'/')
        with open('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'.fa','w+') as f:
            for spe in species_list:
                f.write('>'+spe+paralog1+'\n')
                f.write(Pair_to_sequence[pair][paralog1][spe]+'\n')
                f.write('>'+spe+paralog2+'\n')
                f.write(Pair_to_sequence[pair][paralog2][spe]+'\n')
            f.write('>'+'kluyveri'+paralog1+'\n')
            f.write(Pair_to_sequence[pair]['kluyveri']+'\n')


    
###
#Now output pairs list
###

    with open('/Users/xji3/Genconv/NewClusterPackRun/NewPairsAlignment/pairs_list.txt', 'w+') as f:
        for pair in Pairs_for_alignment:
            f.write('_'.join(pair) + '\n')
