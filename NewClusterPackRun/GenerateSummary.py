from CodonGeneconv import *
import os.path

def get_summary(p_file, output_label = False):
    res = pickle.load(open(p_file, 'r'))
    t = CodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'])
    if res['clock']:
        t.update_by_x_clock(res['x_clock'], Force = res['Force'])
    else:
        t.update_by_x(res['x'], Force = res['Force'])

    nEdge = len(t.edge_to_blen)  # number of edges
    l = nEdge / 2 + 1               # number of leaves
    k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

    leaf_branch = [edge for edge in t.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
    out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
    internal_branch = [x for x in t.edge_to_blen.keys() if not x in leaf_branch]
    assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

    leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    
    out = [t.nsites, res['ll']]
    out.extend(t.pi)
    
    if t.Model == 'HKY': # HKY model doesn't have omega parameter
        out.extend([t.kappa, t.tau])
        label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau']
    elif t.Model == 'MG94':
        out.extend([t.kappa, t.omega, t.tau])
        label = ['length', 'll','pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'omega', 'tau']

    k = len(label)  # record the length of non-blen parameters

    for i in range(len(internal_branch)):
        label.extend([internal_branch[i], leaf_branch[i]])
    label.extend(leaf_branch[i + 1:])
    
##    if t.Model == 'HKY':
##        out.extend([t.edge_to_blen[label[j]] for j in range(k, len(label))])
##    elif t.Model == 'MG94':
##        out.extend([t.edge_to_blen[label[j]] for j in range(k, len(label))])
    out.extend([t.edge_to_blen[label[j]] for j in range(k, len(label))])

    for i in range(k, len(label)):
        label[i] = '(' + ','.join(label[i]) + ')'

    if output_label:
        return out, label
    else:
        return out

def get_expectednum_summary(p_file, output_label = False):

def get_mass_summary(pairs, pair_path, model, summary_file, unfinished_list_file, clock = True, force = False):
    summary_mat = []
    unfinished_list = []
    finished_list = []

    if force:
        prefix = pair_path + 'Force_' + model + '_'
    else:
        prefix = pair_path + model + '_'

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
    p_file = './NewPackageNewRun/HKY_YAL056W_YOR371C_clock.p'
    p_file = './NewPackageNewRun/HKY_YLR406C_YDL075W_clock.p'
    p_file = './NewPackageNewRun/Force_MG94_YML026C_YDR450W_clock.p'
    p_file = './NewPackageNewRun/MG94_YML026C_YDR450W_clock.p'
    res = pickle.load(open(p_file, 'r'))
    print res.keys()
    #print 'X = ', res['x']
    t = CodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'])
    if res['clock']:
        t.update_by_x_clock(res['x_clock'], Force = res['Force'])
    else:
        t.update_by_x(res['x'], Force = res['Force'])

    print 'LL = ', t._loglikelihood(), res['ll']

    print get_summary(p_file, True)

    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])
##
####################################################################################################################################################
##
##    # HKY clock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/HKY_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/HKY_clock_unfinished.txt', clock = True, force = False)
##
##    # HKY nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/HKY_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/HKY_nonclock_unfinished.txt', clock = False, force = False)
##
##    # HKY clock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Force_HKY_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_HKY_clock_unfinished.txt', clock = True, force = True)
##
##    # HKY nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'HKY', summary_file = './NewPackageNewRun/Force_HKY_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_HKY_nonclock_unfinished.txt', clock = False, force = True)
##                
####################################################################################################################################################
##
##    # MG94 clock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/MG94_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/MG94_clock_unfinished.txt', clock = True, force = False)
##
##    # MG94 nonclock model 
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/MG94_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/MG94_nonclock_unfinished.txt', clock = False, force = False)
##
##    # MG94 clock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Force_MG94_clock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_MG94_clock_unfinished.txt', clock = True, force = True)
##
##    # MG94 nonclock model (force tau = 0)
##    get_mass_summary(pairs, pair_path = './NewPackageNewRun/', model = 'MG94', summary_file = './NewPackageNewRun/Force_MG94_nonclock_summary.txt',
##                     unfinished_list_file = './NewPackageNewRun/Force_MG94_nonclock_unfinished.txt', clock = False, force = True)
##                
