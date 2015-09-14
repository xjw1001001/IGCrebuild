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



def _method_1(M, P, R):
    # This is a transcription of method 2 but with vectorization removed.
    #
    nstates, nsites = M.shape
    for site in range(nsites):
        J = np.empty((nstates, nstates))
        for i in range(nstates):
            d = np.zeros(nstates)
            for j in range(nstates):
                d[j] = P[i, j] * R[j, site]
            d = d / d.sum()
            for j in range(nstates):
                J[i, j] = M[i, site] * d[j]
        return J


def _method_2(M, P, R):
    # This uses the original notation.
    #
    nstates, nsites = M.shape
    for site in range(nsites):
        J = np.empty((nstates, nstates))
        for i in range(nstates):
            d = P[i] * R[:, site]
            total = d.sum()
            if total:
                d /= total
            else:
                d = np.zeros_like(d)
            p = M[i, site]
            J[i] = p * d
        return J


def _method_3(M, P, R):
    # Attempt to use a more vectorized notation.
    #
    nstates, nsites = M.shape
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        rdiag = np.diag(r)
        X = P.dot(rdiag)
        y = X.dot(np.ones(nstates))
        ydiag = np.diag(y)
        X = inv(ydiag).dot(X)
        J = np.diag(m).dot(X)
        return J


def _method_4(M, P, R):
    # Attempt to use a more vectorized notation.
    #
    nstates, nsites = M.shape
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        # d = diag(m) * inv(diag(P * r)) * P * diag(r)
        J = np.diag(m).dot(inv(np.diag(P.dot(r)))).dot(P).dot(np.diag(r))
        return J


def _method_5(M, P, R):
    # Attempt to use a more vectorized notation.
    #
    # Note that this implies that the posterior joint distribution
    # divided by the prior conditional distribution
    # is a rank-1 outer product of
    # 1) the marginal initial distribution
    # divided by the matrix-vector product of the prior conditional
    # distribution and the conditional remaining data likelihood and
    # 2) the conditional remaining data likelihood
    #
    nstates, nsites = M.shape
    for site in range(nsites):
        m = M[:, site]
        r = R[:, site]

        # d = diag(m) * inv(diag(P * r)) * P * diag(r)
        J = ((P * r).T * (m / P.dot(r))).T
        return J



def test_joint_distribution_notations():
    nstates = 5
    nsites = 3
    M = np.random.randn(nstates, nsites)
    P = np.random.randn(nstates, nstates)
    R = np.random.randn(nstates, nsites)

    J1 = _method_1(M, P, R)
    J2 = _method_2(M, P, R)
    J3 = _method_3(M, P, R)
    J4 = _method_4(M, P, R)
    J5 = _method_5(M, P, R)

    assert_allclose(J1, J2)
    assert_allclose(J2, J3)
    assert_allclose(J3, J4)
    assert_allclose(J4, J5)
