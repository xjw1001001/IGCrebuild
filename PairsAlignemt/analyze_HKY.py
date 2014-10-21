import sys, os
sys.path.append('/Users/xji3/Genconv')
from CodonBased2Repeats import *

if __name__=='__main__':
    #os.chdir('./PairsAlignemt')

    finished_pair_list_dir = './FinishedPairs.txt'
    pair_list=[]
    with open(finished_pair_list_dir,'r') as f:
        for line in f:
            paralog1 = line[0:7]
            paralog2 = line[8:-1]
            pair_list.append([paralog1,paralog2])

    output_dir = './finished_estimates.txt'
    output_blen_dir = './finished_blen_estimates.txt'
    results_Codon_list = ['RootedTest_Codon_correctedYHR106W_YDR353W_Force_Tau.txt','RootedTest_Codon_correctedYHR106W_YDR353W_Free_Tau.txt'
                    'UnRootedTest_Codon_correctedYHR106W_YDR353W_Force_Tau.txt','UnRootedTest_Codon_correctedYHR106W_YDR353W_Free_Tau.txt']

    results_HKY_list = ['RootedTest_HKYYHR106W_YDR353W_Force_Tau.txt','RootedTest_HKYYHR106W_YDR353W_Free_Tau.txt'
                    'UnRootedTest_HKYYHR106W_YDR353W_Force_Tau.txt','UnRootedTest_HKYYHR106W_YDR353W_Free_Tau.txt']
    
    with open(output_dir,'w+') as f:
        with open(output_blen_dir, 'w+') as ff:
            
            first_line = '\t'.join(['pair','length',
                                  'R0_logL','RF_logL','U0_logL','UF_logL',
                                  'R0_kappa','RF_kappa','U0_kappa','UF_kappa',
                                  'R_Tau','U_Tau',
                                  'R0_TL','RF_TL','U0_TL','UF_TL'])
            first_line_ff = '\t'.join(['pair','length',
                                    '''('N0', 'N1')''', '''('N5', 'cerevisiae')''',
                                    '''('N2', 'N3')''', '''('root', 'kluyveri')''', '''('N4', 'N5')''',
                                    '''('N3', 'kudriavzevii')''', '''('N2', 'bayanus')''', '''('N1', 'N2')''',
                                    '''('N1', 'castellii')''', '''('N5', 'paradoxus')''','''('N4', 'mikatae')''',
                                    '''('N3', 'N4')''', '''('root', 'N0')''','''('N0', 'kluyveri')'''
                                    ])
            ff.write( first_line_ff + '\n' )
            f.write(first_line+'\n')
            for pair in pair_list:
                paralog1 = pair[0]
                paralog2 = pair[1]
                paralog_pair = pair
                numLeaf = 7
                tree_newick = './YeastTree.newick'
                dataloc = './' + paralog1 + '_' + paralog2 + '/' + paralog1 + '_' + paralog2 + '_input.fasta'
                Model = 'Codon_corrected'
                Model = 'HKY'

                with open('./'+paralog1+'_'+paralog2+'/'+paralog1+'_'+paralog2+'_input.fasta','r') as g:
                    g.readline()
                    length = len(g.readline().replace('\n',''))

                Rooted_force_dir = './'+paralog1+'_'+paralog2+'/RootedTest_'+ Model +paralog1+'_'+paralog2+'_Force_Tau.txt'
                with open(Rooted_force_dir,'r') as g:  # Rooted Force Tau output
                    blen = np.ones([2*numLeaf-2])*2
                    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri') ,
                                            cdmodel = True, removegaps = True)
                    parser  = g.readlines()
                    R0_modelnum = int(parser[0].split()[-1])
                    R0_Tau  = float(parser[2].split()[-1])
                    R0_kappa= float(parser[1].split()[-1])
                    R0_blen = [float(parser[c].split()[-1]) for c in range(3,16)]
                    R0_TL   = sum(R0_blen)
                    R0_logL = float(parser[19].split()[-1])
                    R0_para = [float(parser[1].split()[c]) for c in range(3,len(parser[1].split()))]

                    test.Tao = R0_Tau
                    test.para = R0_para
                    test.modelnum = R0_modelnum
                    test.root = 'root'
                    new_root_branch = ('root', test.root_branch[1])
                    test.treetopo.add_node('root')
                    test.treetopo.remove_edge(*test.root_branch)
                    test.treetopo.add_edges_from([('root',test.root_branch[0]),('root',test.root_branch[1])])
                    test.edge_to_blen.pop(test.root_branch)
                    test.edge_to_blen[('root',test.root_branch[0])] = test.blen[0]
                    test.edge_to_blen[('root',test.root_branch[1])] = test.blen[0]
                    test.root_branch = new_root_branch
                    test.blen = R0_blen
                    for c in range(len(test.edge_to_blen)):
                        test.edge_to_blen[test.edge_to_blen.keys()[c]] = test.blen[c]

                    test.save_to_file(Rooted_force_dir.replace('.txt','_class.p'))

                    

                Rooted_free_dir = './'+paralog1+'_'+paralog2+'/RootedTest_'+ Model +paralog1+'_'+paralog2+'_Free_Tau.txt'
                with open(Rooted_free_dir,'r') as g:  # Rooted Free Tau output
                    blen = np.ones([2*numLeaf-2])*2
                    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri'), 
                                            cdmodel = True, removegaps = True)
                    parser = g.readlines()
                    RF_modelnum = int(parser[0].split()[-1])
                    RF_Tau = float(parser[2].split()[-1])
                    RF_kappa= float(parser[1].split()[-1])
                    RF_blen= [float(parser[c].split()[-1]) for c in range(3,16)]
                    RF_TL  = sum(RF_blen)
                    RF_logL= float(parser[19].split()[-1])
                    RF_para = [float(parser[1].split()[c]) for c in range(3,len(parser[1].split()))]

                    test.Tao = RF_Tau
                    test.para = RF_para
                    test.modelnum = RF_modelnum
                    test.root = 'root'
                    new_root_branch = ('root', test.root_branch[1])
                    test.treetopo.add_node('root')
                    test.treetopo.remove_edge(*test.root_branch)
                    test.treetopo.add_edges_from([('root',test.root_branch[0]),('root',test.root_branch[1])])
                    test.edge_to_blen.pop(test.root_branch)
                    test.edge_to_blen[('root',test.root_branch[0])] = test.blen[0]
                    test.edge_to_blen[('root',test.root_branch[1])] = test.blen[0]
                    test.root_branch = new_root_branch
                    test.blen = RF_blen
                    for c in range(len(test.edge_to_blen)):
                        test.edge_to_blen[test.edge_to_blen.keys()[c]] = test.blen[c]

                    test.save_to_file(Rooted_free_dir.replace('.txt','_class.p'))
                               
                UnRooted_force_dir = './'+paralog1+'_'+paralog2+'/UnRootedTest_'+ Model +paralog1+'_'+paralog2+'_Force_Tau.txt'
                with open(UnRooted_force_dir,'r') as g:  # UnRooted Force Tau output
                    blen = np.ones([2*numLeaf-2])*2
                    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri'),
                                            cdmodel = True, removegaps = True)
                    parser = g.readlines()
                    U0_modelnum = int(parser[0].split()[-1])
                    U0_Tau = float(parser[2].split()[-1])
                    U0_kappa= float(parser[1].split()[-1])
                    U0_blen= [float(parser[c].split()[-1]) for c in range(3,15)]
                    U0_TL  = sum(U0_blen)
                    U0_logL= float(parser[18].split()[-1])
                    U0_para = [float(parser[1].split()[c]) for c in range(3,len(parser[1].split()))]

                    test.Tao = U0_Tau
                    test.para = U0_para
                    test.modelnum = U0_modelnum
                    test.blen = U0_blen
                    for c in range(len(test.edge_to_blen)):
                        test.edge_to_blen[test.edge_to_blen.keys()[c]] = test.blen[c]
                    test.save_to_file(UnRooted_force_dir.replace('.txt','_class.p'))
                    

                UnRooted_free_dir = './'+paralog1+'_'+paralog2+'/UnRootedTest_'+ Model +paralog1+'_'+paralog2+'_Free_Tau.txt'
                with open(UnRooted_free_dir,'r') as g:  # Rooted Force Tau output
                    blen = np.ones([2*numLeaf-2])*2
                    test = Codon2RepeatsPhy(numLeaf,blen,tree_newick,dataloc, paralog = paralog_pair, root_branch = ('N0','kluyveri'),
                                            cdmodel = True, removegaps = True)
                    parser = g.readlines()
                    UF_modelnum = int(parser[0].split()[-1])
                    UF_Tau = float(parser[2].split()[-1])
                    UF_kappa= float(parser[1].split()[-1])
                    UF_blen= [float(parser[c].split()[-1]) for c in range(3,15)]
                    UF_TL  = sum(UF_blen)
                    UF_logL= float(parser[18].split()[-1])
                    UF_para = [float(parser[1].split()[c]) for c in range(3,len(parser[1].split()))]

                    test.Tao = UF_Tau
                    test.para = UF_para
                    test.modelnum = UF_modelnum
                    test.blen = UF_blen
                    for c in range(len(test.edge_to_blen)):
                        test.edge_to_blen[test.edge_to_blen.keys()[c]] = test.blen[c]
                    test.save_to_file(UnRooted_free_dir.replace('.txt','_class.p'))

                write_line = '\t'.join([paralog1+'_'+paralog2,str(length),
                                        str(R0_logL),str(RF_logL),str(U0_logL),str(UF_logL),
                                        str(R0_kappa),str(RF_kappa),str(U0_kappa),str(UF_kappa),
                                        str(RF_Tau),str(UF_Tau),
                                        str(R0_TL),str(RF_TL),str(U0_TL),str(UF_TL)])
                f.write(write_line+'\n')
