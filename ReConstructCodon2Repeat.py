import pickle
from CodonBased2Repeats import *
import os

if __name__ == '__main__':
    #os.chdir('/Users/xji3/PairsAlignment/PairsAlignemt')
    finished_pair_list_dir = './FinishedPairs.txt'
    pair_list=[]
    with open(finished_pair_list_dir,'r') as f:
        for line in f:
            paralog1 = line[0:7]
            paralog2 = line[8:-1]
            pair_list.append([paralog1,paralog2])


    for pair in pair_list:
        paralog1 = pair[0]
        paralog2 = pair[1]
        paralog_pair = pair
        ######################################### Rooted Force Tau Case #################################################
        Rooted_force_dir = './'+paralog1+'_'+paralog2+\
                           '/RootedTest_Codon_corrected'+paralog1+'_'+paralog2+'_Force_Tau_class.p'
        load_R0 = pickle.load(open(Rooted_force_dir,'rb'))
        R0_Egeneconv, R0_Enongenconv = load_R0.get_edge_to_expectednumchanges()
        pickle.dump(R0_Egeneconv, open(Rooted_force_dir.replace('_class.p','_Expected_Geneconv_Numbers.p')))
        pickle.dump(R0_Enongeneconv, open(Rooted_force_dir.replace('_class.p','_Expected_NoneGeneconv_Numbers.p')))


        ######################################### Rooted Free Tau Case #################################################
        Rooted_free_dir = './'+paralog1+'_'+paralog2+\
                          '/RootedTest_Codon_corrected'+paralog1+'_'+paralog2+'_Free_Tau_class.p'
        load_Rf = pickle.load(open(Rooted_free_dir,'rb'))
        Rf_Egeneconv, Rf_Enongenconv = load_Rf.get_edge_to_expectednumchanges()
        pickle.dump(Rf_Egeneconv, open(Rooted_free_dir.replace('_class.p','_Expected_Geneconv_Numbers.p')))
        pickle.dump(Rf_Enongeneconv, open(Rooted_free_dir.replace('_class.p','_Expected_NoneGeneconv_Numbers.p')))


        ######################################### UnRooted Force Tau Case #################################################
        UnRooted_force_dir = './'+paralog1+'_'+paralog2+\
                             '/UnRootedTest_Codon_corrected'+paralog1+'_'+paralog2+'_Force_Tau_class.p'
        load_U0 = pickle.load(open(UnRooted_force_dir,'rb'))
        U0_Egeneconv, U0_Enongenconv = load_U0.get_edge_to_expectednumchanges()
        pickle.dump(U0_Egeneconv, open(UnRooted_force_dir.replace('_class.p','_Expected_Geneconv_Numbers.p')))
        pickle.dump(U0_Enongeneconv, open(UnRooted_force_dir.replace('_class.p','_Expected_NoneGeneconv_Numbers.p')))

        
        ######################################### UnRooted Free Tau Case #################################################
        UnRooted_free_dir = './'+paralog1+'_'+paralog2+\
                            '/UnRootedTest_Codon_corrected'+paralog1+'_'+paralog2+'_Free_Tau_class.p'
        load_Uf = pickle.load(open(UnRooted_free_dir,'rb'))
        Uf_Egeneconv, Uf_Enongenconv = load_Uf.get_edge_to_expectednumchanges()
        pickle.dump(Uf_Egeneconv, open(UnRooted_free_dir.replace('_class.p','_Expected_Geneconv_Numbers.p')))
        pickle.dump(Uf_Enongeneconv, open(UnRooted_free_dir.replace('_class.p','_Expected_NoneGeneconv_Numbers.p')))
