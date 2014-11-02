"""
Check equivalent notations for marginal distribution calculations.

M : (nstates, nsites), in
P : (nstates, nstates), in
R : (nstates, nsites), in
D : (nstates, nsites), out

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.linalg import inv


def _method_1(M, P, R):
    # This is a transcription of method 2 but with vectorization removed.
    #
    nstates, nsites = M.shape
    D = np.zeros_like(M)
    for site in range(nsites):
        for i in range(nstates):
            d = np.zeros(nstates)
            for j in range(nstates):
                d[j] = P[i, j] * R[j, site]
            d = d / d.sum()
            for j in range(nstates):
                D[j, site] += M[i, site] * d[j]
    return D


def _method_2(M, P, R):
    # This method uses some vectorization but is not fully vectorized.
    # It is the original implementation.
    #
    nstates, nsites = M.shape
    D = np.zeros_like(M)
    for site in range(nsites):
        for i in range(nstates):
            d = P[i] * R[:, site]
            total = d.sum()
            if total:
                d /= total
            else:
                d = np.zeros_like(d)
            p = M[i, site]
            D[:, site] += p * d
    return D


def _method_3(M, P, R):
    # This is an attempt to more fully vectorize method 2.
    # In particular, the matrix P should be used only through
    # its matrix products with smaller matrices.
    #
    nstates, nsites = M.shape
    D = np.zeros_like(M)
    for site in range(nsites):

        # In this method, use matrix algebra at the per-site level,
        # without combining multiple sites into the formulas.
        m = M[:, site]
        r = R[:, site]
        d = np.zeros(nstates)
        assert_equal(m.shape, (nstates, ))
        assert_equal(r.shape, (nstates, ))

        # d = m.T * inv(diag(P * diag(r) * e)) * P * diag(r)
        rdiag = np.diag(r)
        X = P.dot(rdiag)
        y = X.dot(np.ones(nstates))
        ydiag = np.diag(y)
        X = inv(ydiag).dot(X)
        d = m.dot(X)

        D[:, site] = d

    return D


def _method_4(M, P, R):
    # This is an attempt to more fully vectorize method 3.
    # In particular, the matrix P should be used only through
    # its matrix products with smaller matrices.
    #
    nstates, nsites = M.shape
    D = np.zeros_like(M)

    for site in range(nsites):

        # In this method, use matrix algebra at the per-site level,
        # without combining multiple sites into the formulas.
        m = M[:, site]
        r = R[:, site]
        d = np.zeros(nstates)
        assert_equal(m.shape, (nstates, ))
        assert_equal(r.shape, (nstates, ))

        # The first equation is a direct transcription of
        # method 3 in matrix notation.
        # The second equivalent equation can be written using
        # P only through its forward and adjoint matrix-vector products.
        # d = m.T * inv(diag(P * diag(r) * e)) * P * diag(r)
        # d = m.T * inv(diag(P * r)) * P * diag(r)
        d = m.dot(inv(np.diag(P.dot(r)))).dot(P).dot(np.diag(r))

        D[:, site] = d

    return D


def _method_5(M, P, R):
    # This is an attempt to more fully vectorize method 4.
    # Try removing the explicit loop over sites.
    #
    nstates, nsites = M.shape
    D = np.zeros_like(M)

    """
    for site in range(nsites):

        # In this method, use matrix algebra at the per-site level,
        # without combining multiple sites into the formulas.
        m = M[:, site]
        r = R[:, site]
        d = np.zeros(nstates)
        assert_equal(m.shape, (nstates, ))
        assert_equal(r.shape, (nstates, ))

        # d = m.T * inv(diag(P * r)) * P * diag(r)
        a = m / P.dot(r)
        b = a.dot(P)
        d = b * r

        D[:, site] = d
    """

    A = M / P.dot(R)
    B = A.T.dot(P)
    D = B.T * R

    return D



def test_marginal_distribution_notations():
    nstates = 5
    nsites = 3
    M = np.random.randn(nstates, nsites)
    P = np.random.randn(nstates, nstates)
    R = np.random.randn(nstates, nsites)

    D1 = _method_1(M, P, R)
    D2 = _method_2(M, P, R)
    D3 = _method_3(M, P, R)
    D4 = _method_4(M, P, R)
    D5 = _method_5(M, P, R)

    assert_allclose(D1, D2)
    assert_allclose(D2, D3)
    assert_allclose(D3, D4)
    assert_allclose(D4, D5)
