from CodonGeneconv import *

if __name__ == '__main__':
    
    alignment_file = '/Users/xji3/Genconv/PairsAlignemt/YAL056W_YOR371C/YAL056W_YOR371C_input.fasta'
    newicktree = '/Users/xji3/Genconv/PairsAlignemt/YeastTree.newick'
    paralog = ['YAL056W', 'YOR371C']
    #test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'HKY')#, nnsites = 10)
    test = CodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94')#, nnsites = 10)

    x = np.array([-0.52856124, -0.50895601, -1.08978484,  0.14326525, -2.45656832,
       -1.62703482,  0.47515258, -1.24235701, -0.62424348, -1.00954971,
       -0.80980144, -0.54295938,  0.86255682,  1.21904872, -0.95233185,
       -0.96061732, -0.88000097,  0.57393683])
    test.update_by_x(x)


    print 'Now calculate likelihood'
    #j_ll = test.loglikelihood_and_gradient()

    result = test.get_mle(display = True)
