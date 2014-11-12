"""
Check certain identifiability properties for marginal processes.

Consider a bivariate continuous-time Markov process
for which only one variable is observable.
Furthermore, perhaps that variable is observable
only at the endpoints of an edge.

"""
from __future__ import division, print_function, absolute_import

from itertools import product, permutations

import numpy as np
from numpy.testing import (assert_allclose, assert_equal,
        assert_array_less, assert_)
import scipy.linalg

from jsonctmctree.ll import process_json_in
from jsonctmctree.sampling import(
        sample_distn,
        sample_time_reversible_rate_matrix,
        sample_time_nonreversible_rate_matrix,
        _brute_force_equilibrium,
        )

# FIXME
# I understand why the 'gene conversion' dependence parameter is not
# identifiable when only the endpoints of a single branch are observable,
# but I don't understand why it is still not identifiable when the
# gene conversion rate depends on the distinction between
# synonymous vs. nonsynonymous substitution.
# This distinction can be thought of as partitioning the univariate states
# and scaling rates between partitions differently than the
# rates within states in the same partition.

def _assert_square(M):
    assert_equal(len(Q.shape), 2)
    assert_equal(Q.shape[0], Q.shape[1])


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
        extra homogenization rate
    scaling_mask : 2d ndarray (n, n)
        scales rates

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


def test_identifiability():
    n = 4
    state_space_shape = (n, n)
    nstates = np.prod(state_space_shape)
    bivariate_states = list(product(range(n), repeat=2))
    bivariate_state_pairs = list(permutations(bivariate_states, 2))
    nrepeats = 10
    #for mask in None, _get_symmetric_scaling_mask(n):
    for mask in (None, ):
        for fn_sample in (
            sample_time_reversible_rate_matrix,
            sample_time_nonreversible_rate_matrix,
            ):
            for repeat_index in range(nrepeats):
                for scale in 0.1, 0.5, 2.0:
                    Q, d = fn_sample(n)
                    P = scipy.linalg.expm(scale * Q)
                    for tau in 0, 0.1, 2.0:

                        # Compute the bivariate process.
                        Q_bivariate, d_bivariate = _geneconvify(Q, d, tau, mask)
                        P_bivariate = scipy.linalg.expm(scale * Q_bivariate)

                        # Check that for each (i, i) initial state,
                        # the marginal transition matrix is equal
                        # to the univariate transition matrix.
                        P_marginal_0 = np.zeros((n, n), dtype=float)
                        P_marginal_1 = np.zeros((n, n), dtype=float)
                        for a in range(n):
                            s_initial = (a, a)
                            i = np.ravel_multi_index(
                                    s_initial, state_space_shape)
                            for s in bivariate_states:
                                u, v = s
                                j = np.ravel_multi_index(s, state_space_shape)
                                P_marginal_0[a, u] += P_bivariate[i, j]
                                P_marginal_1[a, v] += P_bivariate[i, j]

                        assert_allclose(P_marginal_0, P)
                        assert_allclose(P_marginal_1, P)
