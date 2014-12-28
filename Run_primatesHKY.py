from CodonGeneconv import *

if __name__ == '__main__':
    alignment_file = '/Users/xji3/Genconv/data/cleaned_input_data.fasta'
    newicktree = '/Users/xji3/Genconv/data/input_tree.newick'
    paralog = ['EDN', 'ECP']
    test = CodonGeneconv(newicktree, alignment_file, paralog, Model = 'HKY')
    print 'Now calculate Likelihood'
    result = test.get_mle(clock = False, display = True)
