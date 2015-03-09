from Rewrite_CodonGeneconv import *
import os.path

if __name__ == '__main__':
##    pairs = []
##    with open('../All_Pairs.txt', 'r') as f:
##        for line in f.readlines():
##            pairs.append(line.replace('\n','').split('_'))
##    pairs.remove(['YLR028C', 'YMR120C'])
##
##    #pairs = pairs[-15:]
##    #pairs = [['YDR502C', 'YLR180W']]
##    for pair in pairs:
##        alignment_file = '/Users/xji3/Genconv/MafftAlignment/' + '_'.join(pair) + '/'+ '_'.join(pair) + '_input.fasta'
##        newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
##        paralog = pair
##        Force = {4:0.0}
##        x = np.array([-0.71060742, -0.51134973, -0.90114847,  2.32335336,  0.33647224,
##       -5.01794203, -2.41741848, -3.48126512, -3.23389773, -4.57747962,
##       -3.66450502, -4.96295342, -3.60311728, -5.80374553, -3.32299616,
##       -3.96429684, -4.25932289])
##
##        test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', clock = False)
##        test.update_by_x(x)
####        test.update_by_x(np.array([ -0.65336901,  -0.53423744,  -0.50935808,   2.73491105,
####        11.40314036,  -1.33157526,  -1.20574659,  -2.34159643,
####        -2.54081475,  -5.5634157 ,  -1.13025978,  -7.12698902,
####        -1.21812082,  -8.43856613,  -1.01319987,  -4.31555508,  -4.64749579]))
##        # This guess would hang YDR502C_YLR180W for HKY nonclock case
##        print 'Now Calculate likelihood for pair ', pair
####        p_file = './NewPackageNewRun/HKY_' + '_'.join(pair) + '_nonclock.p'
####        if os.path.isfile(p_file):
####            res = pickle.load(open(p_file, 'r'))
####            t = ReCodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'], clock = res['clock'], Force = res['Force'])
####            test.update_by_x(t.x)
####            
##        try:
##            result = test.get_mle(display = False)
##            #test.get_ExpectedNumGeneconv()
##            test.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/')
##        except:
##            print 
##            print 'Failed for ', pair, ' non clock'
##            print
##
##        test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', clock = True)
##        p_file = './NewPackageNewRun/HKY_' + '_'.join(pair) + '_clock.p'
####        if os.path.isfile(p_file):
####            res = pickle.load(open(p_file, 'r'))
####            t = ReCodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'], clock = res['clock'], Force = res['Force'])
####            test2.update_by_x_clock(t.x_clock)
####            
##        try:
##            result = test2.get_mle(display = False)
##            #test2.get_ExpectedNumGeneconv()
##            test2.save_to_file( path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/')
##        except:
##            print
##            print 'Failed for ', pair, ' clock'
##            print
##            
##        test3 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', clock = False, Force = Force)
##        test3.update_by_x(x = test.x)
##        try:
##            result = test3.get_mle(display = False)
##            #test.get_ExpectedNumGeneconv()
##            test3.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')
##        except:
##            print 
##            print 'Failed for ', pair, ' non clock Force'
##            print
##
##        test4 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', clock = True, Force = Force)
##
##        test4.update_by_x_clock(test2.x_clock)
##        try:
##            result = test4.get_mle(display = False)
##            #test.get_ExpectedNumGeneconv()
##            test4.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')
##        except:
##            print
##            print 'Failed for ', pair, ' clock Force'
##            print

###############################################################################
#          For unfinished pairs
###############################################################################


##    pair = ['YDR099W', 'YER177W']
##    alignment_file = '/Users/xji3/Genconv/MafftAlignment/' + '_'.join(pair) + '/'+ '_'.join(pair) + '_input.fasta'
##    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
##    paralog = pair
##    Force = {4:0.0}
##    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = True)
##    test.get_mle(display = True, derivative = True)
##    test.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')

    pair = ['YDR502C', 'YLR180W']
    alignment_file = '/Users/xji3/Genconv/MafftAlignment/' + '_'.join(pair) + '/'+ '_'.join(pair) + '_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    paralog = pair
    Force = {4:0.0}
    display_derivative = False

    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = False)
    x = np.array([-0.71060742, -0.51134973, -0.90114847,  2.32335336,  0.33647224,
       -5.01794203, -2.41741848, -3.48126512, -3.23389773, -4.57747962,
       -3.66450502, -4.96295342, -3.60311728, -5.80374553, -3.32299616,
       -3.96429684, -4.25932289])
    test.update_by_x(x)
    test.get_mle(display = display_derivative, derivative = True)
    #test.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')   

    test2 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = False)
    test2.update_by_x(test.x)
    test2.get_mle(display = display_derivative, derivative = True)
    #test2.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')

    test3 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = None, clock = True)
    #test3.update_by_x_clock(x_clock)
    test3.get_mle(display = display_derivative, derivative = True)
    #test.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_')   

    test4 = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY', Force = Force, clock = True)
    test4.update_by_x_clock(test3.x_clock)
    test4.get_mle(display = display_derivative, derivative = True)
    #test2.save_to_file(path = '/Users/xji3/Genconv/NewClusterPackRun/NewPackageNewRun/Force_') 
