"""
Brute force test of weighted expectations on edges.

"""
from __future__ import division, print_function, absolute_import

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from jsonctmctree.expect import process_json_in


def test_tiny_model():

    # Define a tree that is a path.
    # The model transitions from an initial state to an absorbing state.
    nnodes = 4
    nedges = nnodes - 1
    
    for r in (0.1, 0.2):

        print('*** tiny model with r =', r, '***')
        j_in = dict(
                node_count = nnodes,
                process_count = 1,
                state_space_shape = [2],
                prior_feasible_states = [[0]],
                prior_distribution = [1.0],
                tree = dict(
                    row = np.arange(nedges).tolist(),
                    col = np.arange(1, nedges+1).tolist(),
                    process = np.zeros(nedges, dtype=int).tolist(),
                    rate = np.ones(nedges, dtype=int).tolist(),
                    ),
                processes = [dict(
                    row = [[0]],
                    col = [[1]],
                    rate = [r],
                    expect = [1],
                    )],
                observable_nodes = [nnodes-1],
                observable_axes = [0],
                iid_observations = [[1]])

        j_out = process_json_in(j_in)
        print(j_out)


def test_simple_model():

    # Define a tree that is a path with four nodes and three edges.
    # Only the end state is observable.
    # 0 --(0)-- 1 --(1)-- 2 --(3)-- 3
    # The rate matrix will have only three states,
    # consisting of an undecided state and two decided states.
    # The rate towards the decided state 2 is twice as fast
    # as the rate towards the decided state 1.

    np.set_printoptions(precision=16)

    nnodes = 4
    nedges = nnodes - 1

    print('*** simple model ***')

    j_in = dict(
            node_count = 4,
            process_count = 1,
            state_space_shape = [3],
            prior_feasible_states = [[0]],
            prior_distribution = [1.0],
            tree = dict(
                row = np.arange(nedges).tolist(),
                col = np.arange(1, nedges+1).tolist(),
                process = np.zeros(nedges, dtype=int).tolist(),
                rate = np.ones(nedges, dtype=int).tolist(),
                ),
            processes = [dict(
                row = [[0], [0]],
                col = [[1], [2]],
                rate = [0.1, 0.2],
                expect = [1.0, 1.0])],
            observable_nodes = [nnodes-1],
            observable_axes = [0],
            iid_observations = [
                [0],
                [1],
                [0],
                [2]])

    # Compute marginal distributions at internal nodes analytically.
    nnodes = 4
    nstates = 3
    Q = np.array([
        [-0.3, 0.1, 0.2],
        [ 0.0, 0.0, 0.0],
        [ 0.0, 0.0, 0.0]], dtype=float)
    P = expm(Q)
    Px = expm_multiply(Q, np.identity(nstates))
    assert_allclose(P, Px)
    print(P)
    W = np.zeros((nnodes, nstates))
    for assign in product(range(nstates), repeat=nnodes):
        assign = list(assign)
        # compute a likelihood
        if assign[0] != 0:
            continue
        if assign[-1] != 2:
            continue
        lk = 1
        for i, j in zip(assign[:-1], assign[1:]):
            lk *= P[i, j]
        for node, state in enumerate(assign):
            W[node, state] += lk
    print('brute force marginal state distributions at nodes:')
    for node in range(nnodes):
        print('node:', node)
        print('distn:', W[node] / W[node].sum())
    print()

    for r in 0.1, 0.2:
        print('*** brute force r = ', r, '***')
        # Compute marginal distributions at internal nodes analytically.
        nnodes = 4
        nstates = 2
        Q = np.array([
            [-r, r],
            [ 0, 0]], dtype=float)
        P = expm(Q)
        Px = expm_multiply(Q, np.identity(nstates))
        assert_allclose(P, Px)
        print(P)
        W = np.zeros((nnodes, nstates))
        for assign in product(range(nstates), repeat=nnodes):
            assign = list(assign)
            # compute a likelihood
            if assign[0] != 0:
                continue
            if assign[-1] != 1:
                continue
            lk = 1
            for i, j in zip(assign[:-1], assign[1:]):
                lk *= P[i, j]
            for node, state in enumerate(assign):
                W[node, state] += lk
        print('brute force marginal state distributions at nodes:')
        for node in range(nnodes):
            print('node:', node)
            print('distn:', W[node] / W[node].sum())
        print()

    j_out = process_json_in(j_in)
    print(j_out)
