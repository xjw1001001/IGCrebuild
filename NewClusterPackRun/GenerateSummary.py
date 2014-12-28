from CodonGeneconv import *

def get_summary(p_file, output_label = False):
    res = pickle.load(open(p_file, 'r'))
    t = CodonGeneconv(res['newicktree'], res['alignment_file'], res['paralog'], Model = res['Model'])
    t.update_by_x(res['x'])

    nEdge = len(t.edge_to_blen)  # number of edges
    l = nEdge / 2 + 1               # number of leaves
    k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

    leaf_branch = [edge for edge in t.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
    out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
    internal_branch = [x for x in t.edge_to_blen.keys() if not x in leaf_branch]
    assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

    leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
    
    out = []
    out.extend(t.pi)
    
    if t.Model == 'HKY': # HKY model doesn't have omega parameter
        out.extend([t.kappa, t.tau])
    elif t.Model == 'MG94':
        out.extend([t.kappa, t.omega, t.tau])

    label = ['pi_a', 'pi_c', 'pi_g', 'pi_t', 'kappa', 'tau']

    for i in range(len(internal_branch)):
        label.extend([internal_branch[i], leaf_branch[i]])
    label.extend(leaf_branch[i + 1:])

    out.extend([t.edge_to_blen[label[j]] for j in range(6, len(label))])

    if output_label:
        return out, label
    else:
        return out


if __name__ == '__main__':
    p_file = './NewPackageNewRun/HKY_YAL056W_YOR371C_clock.p'
    print get_summary(p_file, True)

    pairs = []
    with open('../All_Pairs.txt', 'r') as f:
        for line in f.readlines():
            pairs.append(line.replace('\n','').split('_'))
    pairs.remove(['YLR028C', 'YMR120C'])

    for pair in pairs:
        p_file = './NewPackageNewRun/HKY_' + '_'.join(pair) + '_clock.p'
                
