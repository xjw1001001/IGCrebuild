from __future__ import print_function, division

from StringIO import StringIO

import dendropy

def read_newick(fin):
    # use dendropy to read this newick file
    t = dendropy.Tree(stream=fin, schema='newick')
    nodes = list(t.preorder_node_iter())
    id_to_idx = {id(n) : i for i, n in enumerate(nodes)}
    edges = []
    edge_rates = []
    for dendro_edge in t.preorder_edge_iter():
        if dendro_edge.tail_node and dendro_edge.head_node:
            na = id_to_idx[id(dendro_edge.tail_node)]
            nb = id_to_idx[id(dendro_edge.head_node)]
            edges.append((na, nb))
            edge_rates.append(dendro_edge.length)
    name_to_node = {str(n.taxon) : id_to_idx[id(n)] for n in t.leaf_nodes()}
    return edges, edge_rates, name_to_node

def main():
    newick_filename = 'collapsed.tree.newick'
    with open(newick_filename) as fin:
        lines = fin.readlines()
    edges, edge_rates, name_to_node = read_newick(StringIO(lines[-1]))
    print('edges:')
    print(edges)
    print('edge rates:')
    print(edge_rates)
    print('name_to_node:')
    print(name_to_node)


main()
