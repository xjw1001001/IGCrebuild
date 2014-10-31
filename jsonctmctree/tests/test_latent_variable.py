"""
Check posterior expectations involving co-evolution with a latent variable.

"""
from __future__ import division, print_function, absolute_import

from itertools import product, permutations

import numpy as np
from numpy.testing import assert_equal, assert_array_less

from jsonctmctree.expect import process_json_in

from sympy import Matrix


def _gen_info(states, norm_rate, conv_rate):
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
            yield sa, sb, rate, expect / rate


def _show_jordan_blocks(R):
    assert_equal(len(R.shape), 2)
    assert_equal(R.shape[0], R.shape[1])
    n = R.shape[0]
    entries = R.flatten().tolist()
    m = Matrix(n, n, entries)
    P, J = m.jordan_form()
    print('rate matrix:')
    print(m)
    print('jordan blocks:')
    print(J)


def test_latent_variable():
    
    # This model is bivariate.
    # Both variables evolve deterministically along a path,
    # but each variable may also convert to its neighbor state.
    n = 3
    norm_rate = 0.1
    conv_rate = 0.01
    state_space_shape = [n, n]
    nstates = np.prod(state_space_shape)
    states = list(product(range(n), repeat=2))
    rate_info = list(_gen_info(states, norm_rate, conv_rate))

    #_show_jordan_blocks(np.array([
        #[-1, 1, 0],
        #[0, -1, 1],
        #[1, 0, -1]], dtype=int))

    # Check the jordan form using sympy.
    """
    mstates = list(product(range(3), repeat=2))
    mnstates = len(mstates)
    rmap = dict((s, i) for i, s in enumerate(mstates))
    R = np.zeros((mnstates, mnstates), dtype=int)
    for sa, sb, r, x in _gen_info(mstates, 1, 1):
        print(sa, sb, r)
        R[rmap[sa], rmap[sb]] = r
    R -= np.diag(R.sum(axis=1))
    _show_jordan_blocks(R)
    """

    nnodes = 5
    nedges = nnodes - 1
    nodes = range(nnodes)
    tree = dict(
            row = nodes[:-1],
            col = nodes[1:],
            rate = [1]*nedges,
            process = [0]*nedges,
            )

    row, col, rate, expect = zip(*rate_info)
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

    j_out = process_json_in(j_in)
    site = 0
    expectations = j_out['edge_expectations'][site]
    assert_equal(len(expectations), nedges)
    assert_array_less(0, expectations)
    assert_array_less(1, expectations[1])
