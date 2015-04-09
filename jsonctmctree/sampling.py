"""
This is for testing only.

"""
from __future__ import division, print_function, absolute_import

from itertools import product, permutations

import numpy as np
from numpy.testing import assert_allclose, assert_equal
import scipy.linalg


def assert_square_matrix(M):
    assert_equal(len(M.shape), 2)
    assert_equal(M.shape[0], M.shape[1])

def assert_symmetric_matrix(M):
    assert_square_matrix(M)
    assert_allclose(M, M.T)

def sample_symmetric_rates(n):
    X = np.random.randn(n, n)
    S = np.exp(X + X.T)
    np.fill_diagonal(S, 0)
    return S

def sample_distn(n):
    d = np.exp(np.random.randn(n))
    return d / d.sum()

def sample_time_reversible_rate_matrix(nstates):
    d = sample_distn(nstates)
    Q = sample_symmetric_rates(nstates) * d
    Q -= np.diag(Q.sum(axis=1))
    assert_allclose(Q.sum(axis=1), 0, atol=1e-13)
    assert_allclose(d.dot(Q), 0, atol=1e-13)
    P = np.diag(d)
    assert_symmetric_matrix(P.dot(Q))
    return Q, d

def _brute_force_equilibrium(Q):
    w, V = scipy.linalg.eig(Q, left=True, right=False)
    i = np.argmin(np.abs(w))
    d = V[:, i].real
    d /= d.sum()
    assert_allclose(Q.sum(axis=1), 0, atol=1e-13)
    assert_allclose(d.dot(Q), 0, atol=1e-13)
    return d

def sample_time_nonreversible_rate_matrix(nstates):
    Q = np.exp(np.random.randn(nstates, nstates))
    Q -= np.diag(Q.sum(axis=1))
    d = _brute_force_equilibrium(Q)
    return Q, d


def _get_symmetric_scaling_mask(n):
    S = np.random.randn(n, n)
    S = np.exp(S + S.T)
    S -= np.diag(np.diag(S))
    return S


def _gen_univariate_transitions(sa, sb):
    for a, b in zip(sa, sb):
        if a != b:
            yield a, b


def _geneconvify(Q_in, d_in, tau, scaling_mask=None):
    """
    Combine two univariate processes into a bivariate process.

    Parameters
    ----------
    Q_in : 2d ndarray (n, n)
        univariate rate matrix
    d_in : 1d ndarray (n, )
        equilibrium distribution
    tau : float
        Additive homogenization rate.
    scaling_mask : 2d ndarray (n, n)
        Scales homogenization rates after having added tau.

    Returns
    -------
    Q_out : 2d ndarray (n*n, n*n)
        bivariate rate matrix
    d_out : 1d ndarray (n*n, )
        equilibrium distribution

    """
    assert_equal(len(d_in.shape), 1)
    assert_equal(len(Q_in.shape), 2)
    n = d_in.shape[0]
    state_space_shape = (n, n)
    assert_equal(Q_in.shape, state_space_shape)
    nstates = np.prod(state_space_shape)
    Q_out = np.zeros((nstates, nstates), dtype=float)
    d_out = np.zeros((nstates, ), dtype=float)
    bivariate_states = list(product(range(n), repeat=2))
    bivariate_state_pairs = list(permutations(bivariate_states, 2))
    for sa, sb in bivariate_state_pairs:
        rate = 0
        t = list(_gen_univariate_transitions(sa, sb))
        if len(t) == 1:
            a, b = t[0]
            rate += Q_in[a, b]
            if sb[0] == sb[1]:
                rate += tau
            if scaling_mask is not None:
                rate *= scaling_mask[a, b]
        if rate:
            i = np.ravel_multi_index(sa, state_space_shape)
            j = np.ravel_multi_index(sb, state_space_shape)
            Q_out[i, j] = rate
    Q_out -= np.diag(Q_out.sum(axis=1))
    d_out = _brute_force_equilibrium(Q_out)
    # Check that the marginal equilbrium distributions
    # of the bivariate process are equal to the equilibrium distribution
    # of the underlying process.
    m0 = np.zeros(n, dtype=float)
    m1 = np.zeros(n, dtype=float)
    for s in bivariate_states:
        a, b = s
        p = d_out[np.ravel_multi_index(s, state_space_shape)]
        m0[a] += p
        m1[b] += p
    assert_allclose(m0, d_in)
    assert_allclose(m1, d_in)
    return Q_out, d_out
