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
    Returns a scipy sparse rate matrix without entries on the diagonal.
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


def get_random_sparse_square_matrix(n):
    """
    Returns a scipy sparse matrix.
    The matrix is square with floating-point type.
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


def check_operator_equivalence(L, M):
    # Check properties of the operator, its transpose, and its adjoint.
    # This is restricted to real square operators.
    np.set_printoptions(precision=10, linewidth=1000, threshold=10000)
    n = L.shape[1]
    B = np.random.randn(n, 2)
    for f, g in (L, M), (L.T, M.T), (L.H, M.T):
        # Check the action of the operator.
        print('fB and gB:')
        print(f.dot(B))
        print(g.dot(B))
        assert_allclose(f.dot(B), g.dot(B))
        # Check the 1-norm and the infinity-norm.
        assert_allclose(f.one_norm(), np.linalg.norm(g, 1))
        assert_allclose(f.inf_norm(), np.linalg.norm(g, np.inf))


def test_RdOperator():
    # This is an n x n square operator.
    n = 4
    R = get_random_rate_matrix(n)
    d = np.random.randn(n)

    # Define the dense numpy ndarray.
    M = R.A + np.diag(d)

    # Define the linear operator.
    L = RdOperator(R, d)

    # Test properties of the operator, its transpose, and its adjoint.
    check_operator_equivalence(L, M)


def test_RdcOperator():
    # This is a 2n x 2n square operator.
    n = 4
    R = get_random_rate_matrix(n)
    d = np.random.randn(n)
    c = np.random.randn(n)

    # Define the dense numpy ndarray.
    Q = R.A + np.diag(d)
    M = np.bmat([[Q, np.diag(c)], [np.zeros((n, n)), Q]]).A

    # Define the linear operator.
    L_Rd = RdOperator(R, d)
    L = RdcOperator(L_Rd, c)

    # Test properties of the operator, its transpose, and its adjoint.
    check_operator_equivalence(L, M)


def test_RdCOperator():
    # This is a 2n x 2n square operator.
    print('*** RdC ***')
    n = 4
    R = get_random_rate_matrix(n)
    d = np.random.randn(n)
    C = get_random_sparse_square_matrix(n)

    # Define the dense numpy ndarray.
    Q = R.A + np.diag(d)
    M = np.bmat([[Q, C.A], [np.zeros((n, n)), Q]]).A

    # Define the linear operator.
    L_Rd = RdOperator(R, d)
    L = RdCOperator(L_Rd, C)

    # Test properties of the operator, its transpose, and its adjoint.
    check_operator_equivalence(L, M)
