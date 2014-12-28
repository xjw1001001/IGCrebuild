from CodonGeneconv import *

if __name__ == '__main__':
    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])
    for pair in pairs:
        alignment_file = '/Users/xji3/Genconv/PairsAlignemt/' + '_'.join(pair) + '/'+ '_'.join(pair) + '_input.fasta'
        newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
        paralog = pair
        Force = {4:0.0}
        test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY')
        print 'Now Calculate likelihood for pair ', pair
##        try:
##            result = test.get_mle(clock = False, display = False)
##            test.get_ExpectedNumGeneconv()
##            test.save_to_file(clock = False, path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/')
##        except:
##            print 
##            print 'Failed for ', pair, ' non clock'
##            print
##
##        try:
##            result = test.get_mle(clock = True, display = False)
##            test.get_ExpectedNumGeneconv()
##            test.save_to_file(clock = True, path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/')
##        except:
##            print
##            print 'Failed for ', pair, ' clock'
##            print

        try:
            test.update_by_x(Force = Force)
            result = test.get_mle(clock = False, display = False, Force = Force)
            test.get_ExpectedNumGeneconv()
            test.save_to_file(clock = False, path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')
        except:
            print 
            print 'Failed for ', pair, ' non clock Force'
            print

        try:
            test.update_by_x(Force = Force)
            result = test.get_mle(clock = True, display = False, Force = Force)
            test.get_ExpectedNumGeneconv()
            test.save_to_file(clock = True, path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')
        except:
            print
            print 'Failed for ', pair, ' clock Force'
            print



##        result = test.get_mle(clock = False, display = False)
##        test.get_ExpectedNumGeneconv()
##        test.save_to_file(clock = False, path = '/Users/xji3/Genconv/NewPackageNewRun/')
##
##        result = test.get_mle(clock = True, display = False)
##        test.get_ExpectedNumGeneconv()
##        test.save_to_file(clock = True, path = '/Users/xji3/Genconv/NewPackageNewRun/')
