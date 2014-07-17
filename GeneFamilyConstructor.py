from Bio.Blast import NCBIXML
import subprocess,os, pickle
from Bio import SeqIO, AlignIO
#, SearchIO
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
        with open(output_dir+i+'/'+i+'_protein.fas','w+') as f:
            f.write('>'+i+'\n')
            f.write(data[i]+'\n')
        
        blastsearch(output_dir+i+'/'+i+'_protein.fas' ,db,evalue , outfmt , output_dir+i+'/'+i+'_database_'+db+'_blastp'+'.xml')


def clean_DNA(data_dna,data_protein,input_dir='/Users/xji3/blast/YeastOutput/'):
    new_data_DNA = {}
    for i in data_dna.keys():
        with open(input_dir+i+'/'+i+'_dna.fas','w+') as f:
            f.write('>'+i+'\n')
            f.write(data_dna[i]+'\n')
        with open(input_dir+i+'/'+i+'_protein.fas','w+') as f:
            f.write('>'+i+'\n')
            f.write(data_protein[i]+'\n')

        cmdline = ['/Users/xji3/Downloads/exonerate-2.2.0/src/program/exonerate','--model','protein2genome',input_dir+i+'/'+i+'_protein.fas', input_dir+i+'/'+i+'_dna.fas']
        align_result = subprocess.check_output(cmdline).split('\n')
        result_iter = align_result.__iter__()
        for line in result_iter:
            if line.strip() and line.split() and line.split()[0:2] == ['Target', 'range:']:
                dna_start = int(line.split()[2])
                dna_stop = int(line.split()[4])
                break
        new_data_DNA[i]=data_dna[i][dna_start:dna_stop]
    return new_data_DNA

def get_one_to_one_dict(data,input_dir = '/Users/xji3/blast/YeastOutput/', output_dir = '/Users/xji3/blast/YeastOutput/'):

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
        length_seq1 = len(data[list(pair)[0]])+0.0
        length_seq2 = len(data[list(pair)[1]])+0.0
        if abs(length_seq1-length_seq2)/min(length_seq1,length_seq2)<threshold:
            new_pair_list.append(pair)
    
    print str(len(pair_list)-len(new_pair_list))+' out of '+str(len(pair_list))+' pairs removed'
    print str(len(new_pair_list)),'remained'
    return new_pair_list

'''
/Users/xji3/Downloads/pamlX-1.2+paml4.7-win32/paml4.7/src
PAML dir
'''


def filter_DS_CodeML(data_dna, pair_list, output_dir = '/Users/xji3/blast/', ctl_file = '/Users/xji3/blast/pairs/codeml.ctl.tmp', threshold = 1.05):
    print
    print 'Apply Ks value filter'
    new_pair_list = []
    check_pair_list = []
    for pair in pair_list:
        id1 = list(pair)[0]
        id2 = list(pair)[1]

        if not os.path.exists(output_dir+'pairs/'):
            os.mkdir(output_dir+'pairs/')

        if not os.path.exists(output_dir+'pairs/'+id1+'_'+id2+'/'):
            os.mkdir(output_dir+'pairs/'+id1+'_'+id2+'/')

        if not os.path.exists(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.fas'):
            with open(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.fas','w+') as f:
                f.write('>'+id1+'\n')
                f.write(data_dna[id1]+'\n')
                f.write('>'+id2+'\n')
                f.write(data_dna[id2]+'\n')
                            
        cmdline = ['../muscle','-in',output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.fas', '-out', output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.aln', '-clw']
        #cmdline = ['../clustalo','--in',output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.fas', '--out', output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.aln']
        subprocess.check_output(cmdline)

        align = AlignIO.read(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.aln',"clustal")

        #Now, remove gaps
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

        max_length = len(align[0].seq)/3
        
        if not os.path.exists(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'_modified.phylip'):
            with open(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'_modified.phylip','w+') as f:
                f.write('2'+'      '+str(max_length*3)+'\n')
                f.write(id1+'\n')
                f.write(str(align[0].seq)[0:3*max_length]+'\n')
                f.write(id2+'\n')
                f.write(str(align[1].seq)[0:3*max_length]+'\n')

        if not os.path.exists(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.trees'):
            with open(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.trees','w+') as f:
                f.write('2      1\n')
                f.write('('+id1+','+id2+');')

        if not os.path.exists(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.ctl.tmp'):
            with open(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.ctl.tmp','w+') as f:
                f.write('seqfile = '+output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'_modified.phylip\n')
                f.write('treefile = '+output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.trees\n')
                f.write('\n')
                f.write('outfile = '+output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'_CodonResult\n')

                f.write(open(ctl_file,'r').read())

        cmdline = ['/Users/xji3/Downloads/pamlX-1.2+paml4.7-win32/paml4.7/src/codeml',output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'.ctl.tmp']
        trial=0
        while(trial<1):
            try:
                subprocess.check_output(cmdline,cwd=output_dir+'pairs/'+id1+'_'+id2+'/')
                trial = 50
            except:
                trial +=1
##        cmdline = ['rm','-r','/Users/xji3/Genconv/GeneFamilyConstructor/lnf','/Users/xji3/Genconv/GeneFamilyConstructor/rst','/Users/xji3/Genconv/GeneFamilyConstructor/rst1','/Users/xji3/Genconv/GeneFamilyConstructor/rub']
##        subprocess.check_output(cmdline).split('\n')


        try:
            parse_paml = codeml.read(output_dir+'pairs/'+id1+'_'+id2+'/'+id1+'_'+id2+'_CodonResult')
            dS = parse_paml['NSsites'][0]['parameters']['dS']
            if dS<threshold:
                new_pair_list.append(pair)
        except:
            print '===============Warning==============='
            print 'Please check pair: ', pair
            check_pair_list.append(pair)


            
    print str(len(pair_list)-len(new_pair_list))+' out of '+str(len(pair_list))+' pairs removed'
    print str(len(new_pair_list)),'remained'
    return new_pair_list,check_pair_list

##def ReadingFrame(seq_dna,seq_protein):
##    translated_seq_1 = 

if __name__ == '__main__':        

    input_file_protein = '/Users/xji3/Genconv/data/Yeast/Saccharomyces_cerevisiae_S288C/cerevisiae_orf_trans.fasta'
    input_file_dna = '/Users/xji3/Genconv/data/Yeast/Saccharomyces_cerevisiae_S288C/cerevisiae_orf_genomic.fasta'

    data_protein = get_data(input_file_protein)
    #run_blast_all(data_protein)
    raw_2copyfamily = get_one_to_one_dict(data_protein)
    pair_list = get_pair_list(raw_2copyfamily)
    data_dna = get_data(input_file_dna)

    #data_dna = clean_DNA(data_dna,data_protein)

    not_match = [i for i in data_dna.keys() if not len(data_dna[i])==3*len(data_protein[i]) ]

    pair_list = filter_DNAlength(data_dna,pair_list)
    pair_list, check_list = filter_DS_CodeML(data_dna,pair_list)
    ##trial = 0
    ##while(len(check_list)>0 and trial<10):
    ##    res_pairs, check_list = filter_DS_CodeML(data_dna,check_list)
    ##    pair_list.extend(res_pairs)
    ##    trial += 1

    ##with open('./current_list.txt','wb+') as f:
    ##    pickle.dump(pair_list,f)

    ##with open('./current_list.txt','rb') as f:
    ##    pair_list = pickle.load(f)

    #run_blast_all(data)
