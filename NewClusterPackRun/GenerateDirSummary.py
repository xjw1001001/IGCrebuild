from DirGeneconv import *
import os.path

def get_summary(p_file, output_label = False):
    res = pickle.load(open(p_file, 'r'))
    t = DirGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'], clock = res['clock'], Force = res['Force'])
    if res['clock']:
        t.update_by_x_clock(res['x_clock'])
    else:
        t.update_by_x(res['x'])

   
    out = [t.nsites, res['ll']]
    out.extend(t.pi)
    
    if t.Model == 'HKY': # HKY model doesn't have omega parameter
        out.extend([t.kappa])
        label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau12', 'tau21']
        out.extend(t.tau)
    elif t.Model == 'MG94':
        out.extend([t.kappa, t.omega])
        label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau12', 'tau21']
        out.extend(t.tau)

    k = len(label)  # record the length of non-blen parameters

    label.extend(res['edge_list'])

    out.extend([t.edge_to_blen[label[j]] for j in range(k, len(label))])

    if res.has_key('ExpectedGeneconv'):
        t.ExpectedGeneconv = res['ExpectedGeneconv']
    if res.has_key('ExpectedDwellTime'):
        t.ExpectedDwellTime = res['ExpectedDwellTime']

    if not t.ExpectedGeneconv:
        t.get_ExpectedNumGeneconv()

    if not t.ExpectedDwellTime:
        t.get_ExpectedHetDwellTime()

    label.extend([ (a, b, 'tau') for (a, b) in res['edge_list']])
    out.extend([t.ExpectedGeneconv[i] / (t.edge_to_blen[i] * t.ExpectedDwellTime[i]) if t.ExpectedDwellTime[i] != 0 else 0 for i in res['edge_list']])


    # Now add directional # of geneconv events
    ExpectedDirectionalNumGeneconv = t._ExpectedDirectionalNumGeneconv()
    label.extend([ (a, b, '1->2') for (a, b) in res['edge_list']])
    out.extend([ExpectedDirectionalNumGeneconv[i][0] for i in res['edge_list']])
    label.extend([ (a, b, '2->1') for (a, b) in res['edge_list']])
    out.extend([ExpectedDirectionalNumGeneconv[i][1] for i in res['edge_list']])


    for i in range(k, len(label)):
        label[i] = '(' + ','.join(label[i]) + ')'

    if output_label:
        return out, label
    else:
        return out


def get_mass_summary(pairs, pair_path, model, summary_file, unfinished_list_file, clock = True, force = False):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if force:
        prefix = pair_path + 'Force_Dir_' + model + '_'
    else:
        prefix = pair_path + 'Dir_' + model + '_'

    if clock:
        suffix = '_clock.p'
    else:
        suffix = '_nonclock.p'

    label = []
    for pair in pairs:
        p_file = prefix + '_'.join(pair) + suffix
        if os.path.isfile(p_file):
            res = get_summary(p_file, True)
            summary_mat.append(res[0])
            label = res[1]
            finished_list.append(pair)
        else:
            unfinished_list.append(pair)
            
    with open(unfinished_list_file, 'w+') as g:
        for pair in unfinished_list:
            g.write('_'.join(pair) + '\n')

    t = np.matrix(summary_mat)
    header = ' '.join(['_'.join(pair) for pair in finished_list])  # column labels
    footer = ' '.join(label)  # row labels
    np.savetxt(open(summary_file, 'w+'), t.T, delimiter = ' ', header = header, footer = footer)
        

if __name__ == '__main__':
##    p_file = './NewPackageNewRun/HKY_YAL056W_YOR371C_clock.p'
##    p_file = './NewPackageNewRun/HKY_YLR406C_YDL075W_clock.p'
##    p_file = './NewPackageNewRun/Force_MG94_YML026C_YDR450W_clock.p'
##    p_file = './NewPackageNewRun/MG94_YML026C_YDR450W_clock.p'
##    res = pickle.load(open(p_file, 'r'))
##    print res.keys()
##    #print 'X = ', res['x']
    

    #pairs = [['YDR502C', 'YLR180W']]
    pairs = [['ECP','EDN']]

#          Primate Dataset

################################################################################################################################################
##
##    # HKY clock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'HKY', summary_file = './NewPackageNewRun/HKY_clock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/HKY_clock_unfinished_Primate.txt', clock = True, force = False)
##
##    # HKY nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'HKY', summary_file = './NewPackageNewRun/HKY_nonclock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/HKY_nonclock_unfinished_Primate.txt', clock = False, force = False)
##
##    # HKY clock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'HKY', summary_file = './NewPackageNewRun/Force_HKY_clock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_HKY_clock_unfinished_Primate.txt', clock = True, force = True)
##
##    # HKY nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'HKY', summary_file = './NewPackageNewRun/Force_HKY_nonclock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_HKY_nonclock_unfinished_Primate.txt', clock = False, force = True)
##                
##################################################################################################################################################
##
##    # MG94 clock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'MG94', summary_file = './NewPackageNewRun/MG94_clock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/MG94_clock_unfinished_Primate.txt', clock = True, force = False)
##
##    # MG94 nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'MG94', summary_file = './NewPackageNewRun/MG94_nonclock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/MG94_nonclock_unfinished_Primate.txt', clock = False, force = False)
##
##    # MG94 clock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'MG94', summary_file = './NewPackageNewRun/Force_MG94_clock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_MG94_clock_unfinished_Primate.txt', clock = True, force = True)
##
##    # MG94 nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/OldResults01222015/', model = 'MG94', summary_file = './NewPackageNewRun/Force_MG94_nonclock_summary_Primate.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_MG94_nonclock_unfinished_Primate.txt', clock = False, force = True)
                


#
#          Yeast Datasets
#

    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])

##################################################################################################################################################

##    # HKY clock model
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Dir_HKY_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_HKY_clock_unfinished.txt', clock = True, force = False)
##
##    # HKY nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Dir_HKY_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_HKY_nonclock_unfinished.txt', clock = False, force = False)
##
##    # HKY clock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Dir_Force_HKY_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_Force_HKY_clock_unfinished.txt', clock = True, force = True)
##
##    # HKY nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Dir_Force_HKY_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_Force_HKY_nonclock_unfinished.txt', clock = False, force = True)
                
##################################################################################################################################################

##    # MG94 clock model
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Dir_MG94_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_MG94_clock_unfinished.txt', clock = True, force = False)

##    # MG94 nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Dir_MG94_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_MG94_nonclock_unfinished.txt', clock = False, force = False)
##
    # MG94 clock model (force tau = 0)
    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Dir_Force_MG94_clock_summary.txt',
                     unfinished_list_file = './NewPackageNewRun/Dir_Force_MG94_clock_unfinished.txt', clock = True, force = True)
##
##    # MG94 nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Dir_Force_MG94_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Dir_Force_MG94_nonclock_unfinished.txt', clock = False, force = True)
                

