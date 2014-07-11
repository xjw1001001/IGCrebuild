from Bio.Blast import NCBIXML
import subprocess,os
from Bio import SeqIO
from Bio.Phylo.PAML import codeml

def blastsearch(query ,db,evalue , outfmt , out):
    return subprocess.check_output([
        '/usr/local/ncbi/blast/bin/blastp',
        '-out',out,
        '-outfmt', str(outfmt),
        '-query',query,
        '-db', db,
        '-evalue',str(evalue)]).split('\n')

def get_data(input_file):
    handle = open(input_file,'rU')
    data = {}
    for record in SeqIO.parse(handle,'fasta'):
        data[record.id]=str(record.seq)
    return data

def run_blast_all(data,outfmt = 5, evalue=1e-18,db='cerevisiae_protein',output_dir='/Users/xji3/blast/YeastOutput/'):
    for i in data.keys():
        try:
            os.mkdir(output_dir+i+'/')
        except:
            print
        with open(output_dir+i+'/'+i+'.fas','w+') as f:
            f.write('>'+i+'\n')
            f.write(data[i]+'\n')
        
        blastsearch(output_dir+i+'/'+i+'.fas' ,db,evalue , outfmt , output_dir+i+'/'+i+'_database_'+db+'_blastp'+'.xml')

def get_one_to_one_dict(data,input_dir = '/Users/xji3/blast/YeastOutput/', output_dir = '/Users/xji3/blast/YeastOutput/'):
    try:
        os.mkdir(output_dir+'2copyFamilies/')
    except:
        print
    Twocopyfamily = {}
    for seq_id in data.keys():
        record = NCBIXML.read(open(input_dir+seq_id+'/'+seq_id+'_database_cerevisiae_protein_blastp.xml'))
        if len(record.alignments)==2:
            Twocopyfamily[seq_id] = [str(hit.hit_id) for hit in record.alignments]
    return Twocopyfamily

def get_pair_list(raw_family):
    pair_list = []
    for gene in raw_2copyfamily:
        if (not set(raw_2copyfamily[gene]) in pair_list) and len(raw_2copyfamily[gene])==2:
            pair_list.append(set(raw_2copyfamily[gene]))
    return pair_list

def create_dbs(db_list=['bayanus','castellii','cerevisiae','kluyveri','kudriavzevii','mikatae','paradoxus'],db_dir='/Users/xji3/blast/db/'):
    for spe in db_list:
        subprocess.check_output([
            '/usr/local/ncbi/blast/bin/makeblastdb',
            '-in', db_dir+spe+'_orf_genomic_1000.fasta',
            '-out',spe,
            '-dbtype','prot']
            ).split('\n')

def filter_DNAlength(data,pair_list,threshold=0.2):
    print 'Apply DNA length filter'
    new_pair_list = []
    for pair in pair_list:
        length_seq1 = len(data[list(pair)[0]])
        length_seq2 = len(data[list(pair)[1]])
        if abs(length_seq1-length_seq2)/min(length_seq1,length_seq2)<threshold:
            new_pair_list.append(pair)
    
    print str(len(pair_list)-len(new_pair_list))+' out of '+str(len(pair_list))+' pairs removed'
    print str(len(new_pair_list)),'remained'
    return new_pair_list
            
        
    
input_file_protein = '/Users/xji3/Genconv/data/Yeast/Saccharomyces_cerevisiae_S288C/cerevisiae_orf_trans.fasta'
data_protein = get_data(input_file_protein)
raw_2copyfamily = get_one_to_one_dict(data_protein)
pair_list = get_pair_list(raw_2copyfamily)

input_file_dna = '/Users/xji3/Genconv/data/Yeast/Saccharomyces_cerevisiae_S288C/cerevisiae_orf_genomic_1000.fasta'
data_dna = get_data(input_file_dna)

pair_list = filter_DNAlength(data_dna,pair_list)
#run_blast_all(data)
