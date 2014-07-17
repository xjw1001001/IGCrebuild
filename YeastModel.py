'''
Geneconv code for Yeast Data
'''

from CodonBased2Repeats import *
from Bio.Blast import NCBIXML
import os

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
        '-outfmt', str(outfmt),
        '-query',query,
        '-db', db,
        '-evalue',str(evalue)]).split('\n')

def run_blast(spe_list,pair_list,data,input_dir = '/Users/xji3/blast/YeastOutput/'):
    for spe in spe_list:
        for pair in pair_list:
            for paralog in pair:
                if not os.path.exists(input_dir+paralog+'/'+paralog+'_dna.fas'):
                    with open(input_dir+paralog+'/'+paralog+'_dna.fas','w+') as f:
                        f.write('>'+paralog+'\n')
                        f.write(data_dna['cerevisiae'][paralog]+'\n')
                blastsearch(input_dir+paralog+'/'+paralog+'_dna.fas',spe+'_dna',1e-10,5,input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml')

def construct_input_fasta(data_dna, spe_list, pair, pair_dir = './data/Pairs/'):
    with open(pair_dir+pair[0]+'_'+pair[1]+'/'+pair[0]+'_'+pair[1]+'_input.fasta','w+') as f:
        for spe in spe_list:
            for paralog in pair:
                if data_dna[spe].has_key(paralog):
                    f.write('>'+spe+paralog+'\n')
                    f.write(data_dna[spe][paralog]+'\n')
                else:
                    print spe, 'does not have gene', paralog

def build_ortholog_dict(spe_list,pair_list,input_dir = '/Users/xji3/blast/YeastOutput/'):
    ortholog_dict = {}
    for pair in pair_list:
        for spe in spe_list:
            for paralog in pair:
                record = NCBIXML.read(open(input_dir+paralog+'/'+paralog+'_database_'+spe+'_dna_blastn.xml'))
                if ortholog_dict.has_key(paralog):
                    ortholog_dict[paralog][spe]=[str(hit.hit_id) for hit in record.alignments]
                else:
                    ortholog_dict[paralog]={spe:[str(hit.hit_id) for hit in record.alignments]}

    return ortholog_dict
                

if __name__=='__main__':
    pair_list_dir = './data/Pairs/pairs_from_paper.txt'
    species_list = ['bayanus','castellii','cerevisiae','kluyveri','kudriavzevii','mikatae','paradoxus']
    pair_list = []
    with open(pair_list_dir,'r') as f:
        for line in f:
            pair_list.append(line.split())

    data_dna = {}
    data_protein = {}
    for spe in species_list:
        print 'constructing sequence dict for species:',spe
        data_dna[spe] = get_data('./data/Yeast/Saccharomyces_'+spe+'/'+spe+'_orf_genomic.fasta')
        data_protein[spe] = get_data('./data/Yeast/Saccharomyces_'+spe+'/'+spe+'_orf_trans.fasta')

    run_blast(species_list,pair_list,data_dna)

    ortholog_dict = build_ortholog_dict(species_list,pair_list)

    

##    for pair in pair_list:
##        construct_input_fasta(data_dna, species_list, pair)


    
