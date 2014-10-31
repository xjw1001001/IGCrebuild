"""
Check posterior expectations involving co-evolution with a latent variable.

"""
from __future__ import division, print_function, absolute_import

from itertools import product, permutations

import numpy as np
from numpy.testing import assert_allclose

from jsonctmctree.expect import process_json_in


def test_latent_variable():
    
    # This model is bivariate.
    # Both variables evolve deterministically along a path,
    # but each variable may also convert to its neighbor state.
    n = 3
    norm_rate = 0.1
    conv_rate = 0.01
    state_space_shape = [n, n]
    states = product(range(n), repeat=2)
    rate_info = []
    for trans in permutations(states, 2):
        sa, sb = trans
        ai, aj = sa
        bi, bj = sb
        rate = 0
        expect = 0
        if ai + 1 == bi and aj == bj:
            rate += norm_rate
        if aj + 1 == bj and ai == bi:
            rate += norm_rate
        if ai != aj and bi == ai and bj == ai:
            rate += conv_rate
            expect += conv_rate
        if ai != aj and bi == aj and bj == aj:
            rate += conv_rate
            expect += conv_rate
        if rate:
            rate_info.append((sa, sb, rate, expect / rate))

    nnodes = 5
    nedges = nnodes - 1
    nodes = range(nnodes)
    tree = dict(
            row = nodes[:-1],
            col = nodes[1:],
            rate = [1]*nedges,
            process = [0]*nedges,
            )

    for info in rate_info:
        print(info)

    row, col, rate, expect = zip(*rate_info)
    print(row)
    print(col)
    print(rate)
    print(expect)
    proc = dict(
            row=row,
            col=col,
            rate=rate,
            expect=expect,
            )

    j_in = dict(
            node_count = nnodes,
            process_count = 1,
            state_space_shape = state_space_shape,
            prior_feasible_states = [[0, 0]],
            prior_distribution = [1.0],
            tree=tree,
            processes=[proc],
            observable_nodes = nodes,
            observable_axes = [0]*nnodes,
            iid_observations = [
                [0, 1, 0, 1, 2],
                ]
            )

    j_out = process_json_in(j_in, debug=True)
    print(j_out)
