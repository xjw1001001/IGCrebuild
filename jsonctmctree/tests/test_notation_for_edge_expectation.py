"""
Check equivalent notations for joint distribution calculations.

Ideally we could avoid ever computing this matrix explicitly.

M : (nstates, nsites), in
P : (nstates, nstates), in
R : (nstates, nsites), in
J : (nstates, nstates), out

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.linalg import inv


def _pseudo_reciprocal(A):
    A_filled = np.where(A, A, 1)
    return np.where(A, np.reciprocal(A_filled), 0)


def _method_1(M, P, R, K):
    # A basic method.
    # Ideally we want to avoid ever computing P or K explicitly.
    # This could be difficult...
    #
    nstates, nsites = M.shape
    P_recip = _pseudo_reciprocal(P)
    S = []
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        # d = diag(m) * inv(diag(P * r)) * P * diag(r)
        J = np.diag(m).dot(inv(np.diag(P.dot(r)))).dot(P).dot(np.diag(r))
        s = np.sum(J * K * P_recip)
        S.append(s)
    return S


def _method_2(M, P, R, K):
    # A basic method.
    # Ideally we want to avoid ever computing P or K explicitly.
    # This could be difficult...
    #
    nstates, nsites = M.shape
    P_recip = _pseudo_reciprocal(P)
    S = []
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        lhs = m / P.dot(r)
        rhs = r
        J_div_P = np.outer(lhs, rhs)
        s = np.sum(J_div_P * K)
        S.append(s)

    return S


def _method_3(M, P, R, K):
    # A basic method.
    # Ideally we want to avoid ever computing P or K explicitly.
    # This could be difficult...
    #
    nstates, nsites = M.shape
    S = []
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        lhs = m / P.dot(r)
        rhs = r
        s = lhs.dot(K).dot(rhs)
        S.append(s)

    return S


def _method_4(M, P, R, K):
    # Fully vectorized.
    A = M * _pseudo_reciprocal(P.dot(R))
    return np.diag(A.T.dot(K.dot(R)))


def test_joint_distribution_notations():
    nstates = 5
    nsites = 3
    M = np.random.randn(nstates, nsites)
    P = np.random.randn(nstates, nstates)
    R = np.random.randn(nstates, nsites)
    K = np.random.randn(nstates, nstates)

    S1 = _method_1(M, P, R, K)
    S2 = _method_2(M, P, R, K)
    S3 = _method_3(M, P, R, K)
    S4 = _method_4(M, P, R, K)

    assert_allclose(S1, S2)
    assert_allclose(S2, S3)
    assert_allclose(S3, S4)
