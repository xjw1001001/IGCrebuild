def read_newick(fin):
    """

    Returns
    -------
    T : undirected weighted networkx tree
        Tree with edge weights.
    root_index : integer
        The root node.
    leaf_name_pairs : sequence
        Sequence of (node, name) pairs.

    """
    # use dendropy to read this newick file
    t = dendropy.Tree(stream=fin, schema='newick')
    leaves = t.leaf_nodes()
    nodes = list(t.postorder_node_iter())
    non_leaves = [n for n in nodes if n not in leaves]
    ordered_nodes = leaves + non_leaves
    root_index = len(ordered_nodes) - 1

    # node index lookup
    node_id_to_index = dict((id(n), i) for i, n in enumerate(ordered_nodes))

    # build the networkx tree
    T = nx.Graph()
    edges = list(t.postorder_edge_iter())
    for i, edge in enumerate(edges):
        if edge.head_node and edge.tail_node:
            na = node_id_to_index[id(edge.head_node)]
            nb = node_id_to_index[id(edge.tail_node)]
            T.add_edge(na, nb, weight=edge.length)

    # get a list of (leaf, name) pairs for the table
    leaf_name_pairs = [(i, str(n.taxon)) for i, n in enumerate(leaves)]

    return T, root_index, leaf_name_pairs
