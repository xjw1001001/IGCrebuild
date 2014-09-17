'''
Geneconv code for Yeast Data
'''

from CodonBased2Repeats import *
from Bio.Blast import NCBIXML
import os
from xlrd import open_workbook

def readInitialPairs(xls_file_dir='./data/Pairs/Supplementary_Table_S1.xls'):
    wb = open_workbook(xls_file_dir)
    sheet = wb.sheets()[0]
    col_list = [str(s.value) for s in sheet.row(0)]
    gene_pairs = {}
    for row in range(1,476):
        pair_name = (str(sheet.cell(row,0).value),str(sheet.cell(row,1).value))
        info = {col_list[col]:sheet.cell(row,col).value for col in range(2,sheet.ncols)}
        gene_pairs[pair_name]=info
            
    return gene_pairs

    


def get_data(input_file):
    handle = open(input_file,'rU')
    data = {}
    for record in SeqIO.parse(handle,'fasta'):
        if record.id[0] == 'Y':
            gene_id = record.id
        else:
            for item in record.description.split():
                if item[0]=='Y' and item[-1]==',':
                    gene_id = item[0:-1]
                    break
            gene_id = record.description.split()[0]
        data[gene_id]=str(record.seq)
    return data

def blastsearch(query ,db,evalue , outfmt , out):
    return subprocess.check_output([
        '/usr/local/ncbi/blast/bin/blastn',
        '-out',out,
        '-reward','1',
        '-penalty','-3',
        '-gapopen','5',
        '-gapextend','2',
        '-task','blastn',
        '-outfmt', str(outfmt),
        '-query',query,
        '-db', db,
        '-evalue',str(evalue)]).split('\n')

def run_blast(spe_list,pair_list,data_dna,evalue = 1e-4,enlarged_database = ['bayanus','kudriavzevii','castellii','kluyveri'],input_dir = '/Users/xji3/blast/YeastOutput/'):
    for spe in spe_list:
        for pair in pair_list:
            for paralog in pair:
                if not os.path.exists(input_dir+paralog+'/'):
                    os.mkdir(input_dir+paralog+'/')

                if not paralog in data_dna['cerevisiae']:
                    pair_list.remove(pair)
                    continue
                
                if not os.path.exists(input_dir+paralog+'/'+paralog+'_dna.fas'):
                    with open(input_dir+paralog+'/'+paralog+'_dna.fas','w+') as f:
                        f.write('>'+paralog+'\n')
                        f.write(data_dna['cerevisiae'][paralog]+'\n')

                if not os.path.exists(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml'):
                    if spe in enlarged_database:
                        blastsearch(input_dir+paralog+'/'+paralog+'_dna.fas',spe+'_dna',evalue,5,input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml')
                        blastsearch(input_dir+paralog+'/'+paralog+'_dna.fas',spe+'_dna2',evalue,5,input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn2.xml')
                    else:
                        db_name = spe+'_dna'
                        blastsearch(input_dir+paralog+'/'+paralog+'_dna.fas',db_name,evalue,5,input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml')

def construct_input_fasta(data_dna, spe_list, pair, pair_dir = './data/Pairs/'):
    with open(pair_dir+pair[0]+'_'+pair[1]+'/'+pair[0]+'_'+pair[1]+'_input.fasta','w+') as f:
        for spe in spe_list:
            for paralog in pair:
                if data_dna[spe].has_key(paralog):
                    f.write('>'+spe+paralog+'\n')
                    f.write(data_dna[spe][paralog]+'\n')
                else:
                    print spe, 'does not have gene', paralog

def build_ortholog_dict(spe_list,pair_list,evalue ,input_dir = '/Users/xji3/blast/YeastOutput/'):
    ortholog_dict = {}
    for pair in pair_list:
        for paralog in pair:
            ortholog_dict[paralog]={spe:[[],[]] for spe in spe_list}
            for spe in spe_list:            
                if not os.path.exists(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml'):
                    blastsearch(input_dir+paralog+'/'+paralog+'_dna.fas',spe+'_dna',evalue,5,input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml')
                
                record = NCBIXML.read(open(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml'))
                ortholog_dict[paralog][spe][0].extend([str(hit.hit_id) for hit in record.alignments])
                
                if os.path.exists(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn2.xml'):
                    record2 = NCBIXML.read(open(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn2.xml'))
                    ortholog_dict[paralog][spe][1].extend([str(hit.hit_id) for hit in record2.alignments])
                
    return ortholog_dict

def filter_have_duplicates_5spe_pairs(pair_list,spe_list,ortholog_dict):
    remove_pair_dict = {spe:[] for spe in spe_list}
    for pair in pair_list:
        paralog1 = pair[0]
        paralog2 = pair[1]
        for spe in spe_list:
            if len(ortholog_dict[paralog1][spe][0])+len(ortholog_dict[paralog1][spe][1])==0 or len(ortholog_dict[paralog2][spe][0])+len(ortholog_dict[paralog2][spe][1])==0:
                if pair not in remove_pair_dict[spe]:
                    remove_pair_dict[spe].append(pair)
                continue
            if not (spe =='castellii' or spe == 'kluyveri'):
                tmp_list = deepcopy(ortholog_dict[paralog1][spe][0])
                tmp_list.extend(ortholog_dict[paralog2][spe][0])
                if len(set(tmp_list))<2:
                    tmp_list = deepcopy(ortholog_dict[paralog1][spe][1])
                    tmp_list.extend(ortholog_dict[paralog2][spe][1])
                    if len(set(tmp_list))<2:
                        remove_pair_dict[spe].append(pair)
                        continue
                

    return remove_pair_dict
        
            

        
##
##class YeastDataset:
##    def __init__(self, list_dir):
##        self.listloc = list_dir
##
##
##    def readpairs(self):
##        

if __name__=='__main__':
    pair_list_dir = './data/Pairs/pairs_from_paper.txt'
    species_list = ['bayanus','cerevisiae','kudriavzevii','mikatae','paradoxus','castellii','kluyveri']
    ncbi_taxonomy_id = {
        'castellii':'Naumovozyma_castellii_CBS_4309',
        
        }
    gene_pairs = readInitialPairs()
    pair_list = gene_pairs.keys()
    evalue = 1e-4

    data_dna = {}
    data_protein = {}
    for spe in species_list:
        print 'constructing sequence dict for species:',spe
        data_dna[spe] = get_data('./data/Yeast/Saccharomyces_'+spe+'/'+spe+'_orf_genomic.fasta')
        data_protein[spe] = get_data('./data/Yeast/Saccharomyces_'+spe+'/'+spe+'_orf_trans.fasta')

    pair_list.remove(('YLR310C', 'YLL016W'))
    pair_list.remove(('YDR134C', 'YLR110C'))

    #run_blast(species_list,pair_list,data_dna,evalue = evalue)

    ortholog_dict = build_ortholog_dict(species_list,pair_list, evalue = evalue)

    remove_pair_dict = filter_have_duplicates_5spe_pairs(pair_list,species_list,ortholog_dict)

    remove_pair_list = []

    for spe in species_list:
        print spe, len(remove_pair_dict[spe])
        remove_pair_list.extend(remove_pair_dict[spe])
    remove_pair_list = list(set(remove_pair_list))

    pair_list = [pair for pair in pair_list if pair not in remove_pair_list]


    direct_pair_list = []    
    for pair in pair_list:
        paralog1 = pair[0]
        paralog2 = pair[1]
        add_to_list = True
        for spe in species_list:
            if spe == 'kluyveri':
                continue
            if not (len(ortholog_dict[paralog1][spe])==2 and  len(ortholog_dict[paralog2][spe])==2):
                add_to_list = False
        if add_to_list:
            direct_pair_list.append(pair)

##    for pair in pair_list:
##        construct_input_fasta(data_dna, species_list, pair)


    
