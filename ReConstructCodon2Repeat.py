import pickle
from CodonBased2Repeats import *
import os

if __name__ == '__main__':
    os.chdir('/Users/xji3/PairsAlignment')
    numLeaf = 7
    blen = np.ones([2*numLeaf-2])*2

    tree_newick = './YeastTree.newick'
    dataloc = './YBL087C_YER117W/YBL087C_YER117W_input.fasta'
    paralog_pair = ['YBL087C','YER117W']

    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,cdmodel = True, removegaps = True)

    pickle.dump(test, open('pickle_test.p','wb+'))

    test_load = pickle.load(open('pickle_test.p','rb'))

