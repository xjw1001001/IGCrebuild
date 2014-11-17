"""
Use a laplace transform to solve a puzzle related to the matrix exponential.

"""
from __future__ import division, print_function, absolute_import

from itertools import product, permutations

import numpy as np
from numpy.testing import (assert_allclose, assert_equal,
        assert_array_less, assert_)
import scipy.linalg

from jsonctmctree.util import assert_square
from jsonctmctree.ll import process_json_in
from jsonctmctree.sampling import(
        sample_distn,
        sample_time_reversible_rate_matrix,
        sample_time_nonreversible_rate_matrix,
        _brute_force_equilibrium,
        _geneconvify,
        assert_symmetric_matrix,
        )


def _vec(M):
    return M.T.flatten()


def _identity_like(M):
    assert_equal(len(M.shape), 2)
    assert_equal(M.shape[0], M.shape[1])
    n = M.shape[0]
    return np.identity(n, dtype=M.dtype)

def _kronecker_sum_1d(a, b):
    return _vec(a[np.newaxis, :] + b[:, np.newaxis])

def _kronecker_sum(a, b):
    return np.kron(a, _identity_like(b)) + np.kron(_identity_like(a), b)


def test_kronecker_sum_1d():
    a = np.random.randn(3)
    b = np.random.randn(4)
    actual = _kronecker_sum_1d(a, b)
    desired = np.diag(_kronecker_sum(np.diag(a), np.diag(b)))
    assert_allclose(actual, desired)


def _check_laplace_transform_equilibrium_solution(tau):
    # The tau parameter must be positive for this solution to work.
    assert_array_less(0, tau)

    # Define the structural properties of the process.
    n = 5
    state_space_shape = (n, n)
    nstates = np.prod(state_space_shape)
    bivariate_states = list(product(range(n), repeat=2))
    bivariate_state_pairs = list(permutations(bivariate_states, 2))
    I_univariate = np.identity(n)

    # Sample a random time-reversible univariate rate matrix
    # together with its stationary distribution.
    Q, d = sample_time_reversible_rate_matrix(n)

    # Compute a decomposition of Q.
    a = np.sqrt(d)
    b = np.reciprocal(a)
    S = Q * np.outer(a, b)
    assert_allclose(S, S.T, atol=1e-12)
    w, U = scipy.linalg.eigh(S)
    assert_allclose(U.dot(U.T), I_univariate, atol=1e-12)
    assert_allclose(U.dot(np.diag(w)).dot(U.T), S)

    # Check properties of the Kronecker sum of Q and Q.
    # This corresponds to joint independent evolution of two variables.
    # First check detailed balance.
    # Then check the eigendecomposition.
    R_Q = _kronecker_sum(Q, Q)
    R_d = np.kron(d, d)
    assert_allclose(R_d.sum(), 1)
    assert_symmetric_matrix(np.diag(R_d).dot(R_Q))
    R_U = np.kron(U, U)
    R_w = _kronecker_sum_1d(w, w)
    R_S = R_U.dot(np.diag(R_w)).dot(R_U.T)
    R_a = np.sqrt(R_d)
    R_b = np.reciprocal(R_a)
    assert_allclose(R_S * np.outer(R_b, R_a), R_Q, atol=1e-12)

    # Get the bivariate gene-conversion process rate matrix directly.
    # Define the gene conversion mask.
    T = np.ones((n, n)) - np.identity(n)
    mask = _kronecker_sum(T, T) * _vec(np.identity(n))
    T_Q = R_Q + tau * mask
    T_Q = T_Q - np.diag(T_Q.sum(axis=1))
    T_d_brute = _brute_force_equilibrium(T_Q)
    assert_equal(T_d_brute.shape, (n*n, ))

    # Check the shortcut for computing the equilibrium distribution
    # of the bivariate process with gene conversion.
    # This involves the Laplace transform.
    s = 2 * tau
    T_u = 1 / (1 - R_w / s)
    T_S = R_U.dot(np.diag(T_u)).dot(R_U.T)
    laplace_thing = T_S * np.outer(R_b, R_a)
    T_d_clever = _vec(np.diag(d)).dot(laplace_thing)
    assert_equal(T_d_clever.shape, (n*n, ))
    assert_allclose(T_d_clever, T_d_brute)

    # Check separate calculations of expectations.
    M_u = 1 / (1 - w / tau)
    yet_another_laplace_thing = U.dot(np.diag(M_u)).dot(U.T) * np.outer(b, a)
    T_d_yet_another = np.diag(d).dot(yet_another_laplace_thing).flatten()
    assert_equal(T_d_yet_another.shape, (n*n, ))
    assert_allclose(T_d_yet_another, T_d_brute)

    # Check some transition probabilities.
    P = U.dot(np.diag(np.exp(w))).dot(U.T) * np.outer(b, a)
    R_P = R_U.dot(np.diag(np.exp(R_w))).dot(R_U.T) * np.outer(R_b, R_a)
    desired = np.array([np.kron(r, r) for r in P])
    actual = R_P[_vec(np.identity(n, dtype=bool)), :]
    assert_allclose(actual, desired)



#FIXME under construction
def _check_three_paralogs(tau):
    # The tau parameter must be positive for this solution to work.
    assert_array_less(0, tau)

    # Define the structural properties of the process.
    n = 5
    state_space_shape = (n, n)
    nstates = np.prod(state_space_shape)
    bivariate_states = list(product(range(n), repeat=2))
    bivariate_state_pairs = list(permutations(bivariate_states, 2))
    I_univariate = np.identity(n)

    # Sample a random time-reversible univariate rate matrix
    # together with its stationary distribution.
    Q, d = sample_time_reversible_rate_matrix(n)

    # Compute a decomposition of Q.
    a = np.sqrt(d)
    b = np.reciprocal(a)
    S = Q * np.outer(a, b)
    assert_allclose(S, S.T, atol=1e-12)
    w, U = scipy.linalg.eigh(S)
    assert_allclose(U.dot(U.T), I_univariate, atol=1e-12)
    assert_allclose(U.dot(np.diag(w)).dot(U.T), S)

    # Check properties of the Kronecker sum of Q and Q.
    # This corresponds to joint independent evolution of two variables.
    # First check detailed balance.
    # Then check the eigendecomposition.
    R_Q = _kronecker_sum(Q, Q)
    R_d = np.kron(d, d)
    assert_allclose(R_d.sum(), 1)
    assert_symmetric_matrix(np.diag(R_d).dot(R_Q))
    R_U = np.kron(U, U)
    R_w = _kronecker_sum_1d(w, w)
    R_S = R_U.dot(np.diag(R_w)).dot(R_U.T)
    R_a = np.sqrt(R_d)
    R_b = np.reciprocal(R_a)
    assert_allclose(R_S * np.outer(R_b, R_a), R_Q, atol=1e-12)

    # Get the bivariate gene-conversion process rate matrix directly.
    # Define the gene conversion mask.
    T = np.ones((n, n)) - np.identity(n)
    mask = _kronecker_sum(T, T) * _vec(np.identity(n))
    T_Q = R_Q + tau * mask
    T_Q = T_Q - np.diag(T_Q.sum(axis=1))
    T_d_brute = _brute_force_equilibrium(T_Q)

    # Check the shortcut for computing the equilibrium distribution
    # of the bivariate process with gene conversion.
    # This involves the Laplace transform.
    s = 2 * tau
    T_u = 1 / (1 - R_w / s)
    T_S = R_U.dot(np.diag(T_u)).dot(R_U.T)
    laplace_thing = T_S * np.outer(R_b, R_a)
    T_d_clever = _vec(np.diag(d)).dot(laplace_thing)
    assert_allclose(T_d_clever, T_d_brute)

    # Check separate calculations of expectations.
    M_u = 1 / (1 - w / tau)
    yet_another_laplace_thing = U.dot(np.diag(M_u)).dot(U.T) * np.outer(b, a)
    T_d_yet_another = np.diag(d).dot(yet_another_laplace_thing).flatten()
    assert_allclose(T_d_yet_another, T_d_brute)

    # Redundant test.
    assert_allclose(T_d_clever, T_d_yet_another)

    # Check some transition probabilities.
    P = U.dot(np.diag(np.exp(w))).dot(U.T) * np.outer(b, a)
    R_P = R_U.dot(np.diag(np.exp(R_w))).dot(R_U.T) * np.outer(R_b, R_a)
    desired = np.array([np.kron(r, r) for r in P])
    actual = R_P[_vec(np.identity(n, dtype=bool)), :]
    assert_allclose(actual, desired)


def test_laplace_transform_equilibrium_solution():
    nrepeats = 10
    for repeat_index in range(nrepeats):
        for tau in 0.1, 2.0:
            _check_laplace_transform_equilibrium_solution(tau)


def test_empirical():
    nsamples = 10000
    n = 3
    I_univariate = np.identity(n)

    # Sample a random time-reversible univariate rate matrix
    # together with its stationary distribution.
    Q, d = sample_time_reversible_rate_matrix(n)

    # Compute a decomposition of Q.
    a = np.sqrt(d)
    b = np.reciprocal(a)
    S = Q * np.outer(a, b)
    assert_allclose(S, S.T, atol=1e-12)
    w, U = scipy.linalg.eigh(S)
    assert_allclose(U.dot(U.T), I_univariate, atol=1e-12)
    assert_allclose(U.dot(np.diag(w)).dot(U.T), S)

    # Check properties of the Kronecker sum of Q and Q.
    # This corresponds to joint independent evolution of two variables.
    # First check detailed balance.
    # Then check the eigendecomposition.
    R_Q = _kronecker_sum(Q, Q)
    R_d = np.kron(d, d)
    assert_allclose(R_d.sum(), 1)
    assert_symmetric_matrix(np.diag(R_d).dot(R_Q))
    R_U = np.kron(U, U)
    R_w = _kronecker_sum_1d(w, w)
    R_S = R_U.dot(np.diag(R_w)).dot(R_U.T)
    R_a = np.sqrt(R_d)
    R_b = np.reciprocal(R_a)
    assert_allclose(R_S * np.outer(R_b, R_a), R_Q, atol=1e-12)

    # Sample times independently per branch.
    P_A = np.zeros((n, n*n), dtype=float)
    for sample_index in range(nsamples):
        t0 = np.random.exponential(scale=1)
        P0 = U.dot(np.diag(np.exp(t0*w))).dot(U.T) * np.outer(b, a)
        t1 = np.random.exponential(scale=1)
        P1 = U.dot(np.diag(np.exp(t1*w))).dot(U.T) * np.outer(b, a)
        P_A += np.array([np.kron(P0[i], P1[i]) for i in range(n)])
    P_A /= nsamples

    # Sample one time shared across both branches.
    P_B = np.zeros((n, n*n), dtype=float)
    for sample_index in range(nsamples):
        t0 = np.random.exponential(scale=1)
        P0 = U.dot(np.diag(np.exp(t0*w))).dot(U.T) * np.outer(b, a)
        P_B += np.array([np.kron(P0[i], P0[i]) for i in range(n)])
    P_B /= nsamples

    # Sample one longer time shared across both branches.
    P_C = np.zeros((n, n*n), dtype=float)
    for sample_index in range(nsamples):
        t0 = np.random.exponential(scale=2)
        P0 = U.dot(np.diag(np.exp(t0*w))).dot(U.T) * np.outer(b, a)
        P_C += np.array([np.kron(P0[i], P0[i]) for i in range(n)])
    P_C /= nsamples

    # Sample one shorter time shared across both branches.
    P_D = np.zeros((n, n*n), dtype=float)
    for sample_index in range(nsamples):
        t0 = np.random.exponential(scale=0.5)
        P0 = U.dot(np.diag(np.exp(t0*w))).dot(U.T) * np.outer(b, a)
        P_D += np.array([np.kron(P0[i], P0[i]) for i in range(n)])
    P_D /= nsamples

    # Use the Laplace transform to exactly compute
    # one of the empirically estimated transition matrices.
    P0 = U.dot(np.diag(1 / (1 - w))).dot(U.T) * np.outer(b, a)
    P_E = np.array([np.kron(P0[i], P0[i]) for i in range(n)])
    assert_allclose(P_E.sum(axis=1), 1)

    # Use the Laplace transform to exactly compute
    # one of the empirically estimated transition matrices.
    P0 = R_U.dot(np.diag(1 / (1 - R_w))).dot(R_U.T) * np.outer(R_b, R_a)
    P_F = P0[_vec(np.identity(n, dtype=bool)), :]
    assert_allclose(P_F.sum(axis=1), 1)

    #print('same time scale:')
    #print(P_B - P_A)
    #print(d.dot(P_B) - d.dot(P_A))
    #print()

    #print('double time scale:')
    #print(P_C - P_A)
    #print(d.dot(P_C) - d.dot(P_A))
    #print()

    #print('half time scale:')
    #print(P_D - P_A)
    #print(d.dot(P_D) - d.dot(P_A))
    #print()

    #print('analytic laplace for independent times:')
    #print(P_E - P_A)
    #print()

    #print('analytic laplace for shared times:')
    #print(P_F - P_B)
    #print()
