"""
Given a newick gene tree, return a species tree with averaged branch lengths.

The leaves of the tree should be named according to the taxa.

"""
from __future__ import print_function, division


# Example newick gene tree input.
"""
(((((((cerevisiaeYDR450W:0.0208250416694,paradoxusYDR450W:0.0187907761474):0.00187312520584,mikataeYDR450W:0.0447343695086):0.0101057212007,kudriavzeviiYDR450W:0.0302301795725):0.0135677728727,bayanusYDR450W:0.0372329696241):0.0537401647561,castelliiYDR450W:0.0586859466806):0.0200584641939,(((((cerevisiaeYML026C:0.0208250416694,paradoxusYML026C:0.0187907761474):0.00187312520584,mikataeYML026C:0.0447343695086):0.0101057212007,kudriavzeviiYML026C:0.0302301795725):0.0135677728727,bayanusYML026C:0.0372329696241):0.0537401647561,castelliiYML026C:0.0586859466806):0.0200584641939):0.1,kluyveriYML026C:0.032938469918);
"""


import argparse
import sys

import networkx as nx

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


def process_newick_string(newick_stream_in, newick_stream_out, paralogs):

    # Define the map from paralog name to index.
    paralog_name_to_idx = {s : i for i, s in enumerate(paralogs)}

    # Parse the newick string.
    edges, edge_rates, name_to_node = read_newick(newick_stream_in)

    # Break the full names into species and gene components.
    # At the same time, track the species names,
    # and track each paralog and species index.
    node_to_paralog_idx = dict()
    node_to_species_idx = dict()
    species_name_to_idx = dict()
    species_list = []
    for full_name, node in name_to_node.items():
        paralog_names = []
        species_names = []
        for p in paralogs:
            if full_name.endswith(p):
                paralog_names.append(p)
                species_names.append(full_name[:len(p)])
        if not species_names or not paralog_names:
            raise Exception('The input name "%s" does not seem to be in '
                    'the form of a species name prefix followed by a gene '
                    'name suffix, for the given list of gene names "%s".' % (
                        full_name, paralogs))
        if len(paralog_names) != 1 or len(species_names) != 1:
            raise Exception('Failed to uniquely parse the name "%s" '
                    'into a species name followed by a gene name.' % full_name)
        paralog_name = paralog_names[0]
        species_name = species_names[0]

        # Map the node of the full tree to the species index.
        species_idx = species_name_to_idx.get(species_name, None)
        if species_idx is None:
            species_idx = len(species_list)
            species_name_to_idx[species_name] = species_idx
            species_list.append(species_name)
        node_to_species_idx[node] = species_idx

        # Map the node of the full tree to the paralog index.
        paralog_idx = paralog_name_to_idx.get(paralog_name, None)
        if paralog_idx is None:
            raise Exception('Failed to interpret the paralog name for '
                    'the full name "%s"' % full_name)
        node_to_paralog_idx[node] = paralog_idx

    # Create the full gene tree.
    T = nx.DiGraph()
    T.add_edges_from(edges)

    # For each node of the full gene tree,
    # map the node to a sorted tuple of subtree species indices.
    node_to_subtree_species_indices = dict()
    for na in nx.dfs_postorder_nodes(T):
        subtree_species_indices = set()
        successors = list(T.successors(na))
        if successors:
            for nb in successors:
                subtree_species_indices.update(
                        node_to_subtree_species_indices[nb])
        else:
            subtree_species_indices.add(node_to_species_idx[na])
        node_to_subtree_species_indices[na] = tuple(
                sorted(subtree_species_indices))

    print(len(T))
    print(len(T.edges()))
    print(node_to_subtree_species_indices)
    print('species list:', species_list)


def main(args):
    # Read the newick string defining the gene tree on stdin,
    # and write the newick string defining the species tree on stdout,
    # averaging the corresponding branch lengths.
    process_newick_string(sys.stdin, sys.stdout, args.paralogs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paralogs', nargs='+')
    args = parser.parse_args()
    main(args)
