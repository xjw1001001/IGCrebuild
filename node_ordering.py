"""
Functions related to node ordering.

"""
from __future__ import print_function

import networkx as nx
from numpy.testing import assert_equal


def get_node_to_depth(T, root):
    node_to_depth = {root : 0}
    for head, tail in nx.bfs_edges(T, root):
        node_to_depth[tail] = node_to_depth[head] + 1
    return node_to_depth


def get_node_to_subtree_depth(T, root):
    subdepth = {}
    for node in reversed(list(nx.bfs_nodes(T, root))):
        successors = T.successors(tail)
        if not successors:
            subdepth[node] = 0
        else:
            subdepth[node] = max(subdepth[s] for s in successors) + 1
    return subdepth


def get_node_to_subtree_thickness(T, root):
    thickness = {}
    for node in nx.dfs_postorder_nodes(T, root):
        successors = T.successors(node)
        if not successors:
            thickness[node] = 1
        else:
            w = sorted((thickness[n] for n in successors), reverse=True)
            thickness[node] = max(w + i for i, w in enumerate(w))
    return thickness


def get_node_evaluation_order(T, root):
    thickness = get_node_to_subtree_thickness(T, root)
    expanded = set()
    stack = [root]
    while stack:
        n = stack.pop()
        if n in expanded:
            yield n
        else:
            successors = list(T.successors(n))
            if not successors:
                yield n
            else:
                expanded.add(n)
                pairs = [(thickness[x], x) for x in successors]
                progeny = [x for w, x in sorted(pairs)]
                stack.extend([n] + progeny)


def get_example_tree():
    T = nx.DiGraph()
    T.add_edges_from((
        (0, 1),
        (1, 2),
        (1, 3),
        (0, 4),
        (4, 5),
        (5, 6),
        (5, 7),
        (4, 8)))
    root = 0
    return T, root


def test_subtree_thickness():
    T, root = get_example_tree()
    d_actual = get_node_to_subtree_thickness(T, root)
    d_desired = {
            0 : 3,
            1 : 2,
            2 : 1,
            3 : 1,
            4 : 2,
            5 : 2,
            6 : 1,
            7 : 1,
            8 : 1}
    assert_equal(d_actual, d_desired)


def test_node_evaluation_order():
    T, root = get_example_tree()
    v_actual = list(get_node_evaluation_order(T, root))
    v_desired = (7, 6, 5, 8, 4, 3, 2, 1, 0)
    assert_equal(v_actual, v_desired)


def main():
    test_subtree_thickness()
    test_node_evaluation_order()


if __name__ == '__main__':
    main()
