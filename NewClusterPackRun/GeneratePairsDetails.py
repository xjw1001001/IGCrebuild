from Bio import SeqIO

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
        data[gene_id]=str(record.description)
    return data


if __name__ == '__main__':

    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])

    input_file = '../data/Yeast/Saccharomyces_cerevisiae/cerevisiae_orf_genomic.fasta'
    cerevisiae_data = get_data(input_file)

    with open('./cerevisiae_genes.txt', 'w+') as f:
        with open('./cerevisiae_genes_description.txt', 'w+') as g:
            for pair in pairs:
                for i in range(2):
                    f.write(pair[i] + '\n')
                    g.write(cerevisiae_data[pair[i]] + '\n')
        
