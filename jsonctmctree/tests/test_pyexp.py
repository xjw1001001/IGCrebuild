"""
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal, assert_array_less, assert_allclose

import scipy.sparse

import jsonctmctree
from jsonctmctree.pyexp.ctmc_ops import RdOperator, RdcOperator, RdCOperator

def get_random_rate_matrix(n):
    """
    This is a scipy sparse rate matrix without entries on the diagonal.

    """
    row_ind = []
    col_ind = []
    shape = (n, n)
    mask = np.random.randint(0, 2, size=shape)
    for i in range(n):
        for j in range(n):
            if i != j and mask[i, j]:
                row_ind.append(i)
                col_ind.append(j)
    ndata = len(row_ind)
    data = np.exp(np.random.randn(ndata))
    R = scipy.sparse.csr_matrix((data, (row_ind, col_ind)), shape=shape)
    return R


def test_RdOperator():
    n = 4
    R = get_random_rate_matrix(n)
    d = np.random.randn(n)

    # Define the linear operator.
    L = RdOperator(R, d)

    # Define the dense numpy ndarray.
    M = R.A + np.diag(d)

    # Pick a random matrix for testing.
    B = np.random.randn(n, 2)

    # Check the forward application of the operator.
    assert_allclose(L.dot(B), M.dot(B))

    # Check the action of the adjoint.
    # Because the operator is real, the adjoint is the transpose.
    assert_allclose(L.H.dot(B), M.T.dot(B))
    assert_allclose(L.T.dot(B), M.T.dot(B))
