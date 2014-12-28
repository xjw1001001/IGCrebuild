from CodonGeneconv import *

if __name__ == '__main__':
    pair = ['YLR406C', 'YDL075W']
    alignment_file = '/Users/xji3/Genconv/PairsAlignemt/' + '_'.join(pair) + '/'+ '_'.join(pair) + '_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    paralog = pair
    test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94')

    print 'Now Calculate likelihood for pair ', pair
    clock = False
    result = test.get_mle(clock = clock, display = True)
    test.get_ExpectedNumGeneconv()
    test.save_to_file(clock = clock, path = '/Users/xji3/Genconv/NewPackageNewRun/')

    clock = True
    result = test.get_mle(clock = clock, display = True)
    test.get_ExpectedNumGeneconv()
    test.save_to_file(clock = clock, path = '/Users/xji3/Genconv/NewPackageNewRun/')
