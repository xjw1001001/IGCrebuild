"""
Test functions related to node ordering.

"""
from __future__ import division, print_function, absolute_import

import networkx as nx
from numpy.testing import assert_equal

from jsonctmctree.node_ordering import (
        get_node_to_depth,
        get_node_to_subtree_depth,
        get_node_to_subtree_thickness,
        get_node_evaluation_order,
        )


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
