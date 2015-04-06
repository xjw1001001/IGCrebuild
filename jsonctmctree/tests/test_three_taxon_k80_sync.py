"""
Numerically check identifiability for a somewhat complicated toy system.

This is a 2-taxon system with an out-taxon, an in-taxon and a single gene
duplication with IGC along the duplicated lineage.
The model will be K80 with unknown kappa, so the system will have four
degrees of freedom:
 * the kappa ratio
 * the IGC rate on the post-duplication lineage
 * the out-taxon lineage substitution rate*time
 * the in-taxon lineage substitution rate*time

According to some back of the envelope calculations, various symmetries
reduce the 4^3 = 64 combinations of joint nucleotide states
at aligned sites to only 7 distinct categories.
A distribution over these 7 categories has 6 degrees of freedom,
so a simple counting argument cannot rule out the possibility that
the system is identifiable (identifying 4 process parameters
from 6 pieces of information).

An intuitive argument suggests that it should be identifiable.
The data pair consisting of the out-group and one of the in-group
gene duplicates should be enough to consistently estimate
the kappa ratio and the total substitution rate*time between the out-taxon
and the in-taxon.  The data pair consisting of the two duplicate in-group
genes should be enough to identify the rate*time since duplication
and the synchronization rate, given the previously identified kappa ratio.
So together, this should be able to consistently identify all four of the
parameters in the model.

The multivariate states will be represented by ndarrays.

"""
from __future__ import division, print_function, absolute_import


from functools import partial
from itertools import product

import numpy as np
from numpy.testing import assert_allclose, assert_equal
from scipy.special import kl_div
from scipy.linalg import expm, eigvalsh, inv

KAPPA = 4
TRANSITIONS = set([(0, 1), (1, 0), (2, 3), (3, 2)])


def mravel(indices):
    return np.ravel_multi_index(indices, (4, 4, 4))


def add_rates(Q, kappa, sync_rate, r0, r1):
    for sa in product(range(4), repeat=3):
        sa = np.array(sa)
        idxa = mravel(sa)
        ai, aj, ak = sa

        # Add rates corresponding to substitutions.
        # There are three lineages.
        for variable in range(3):
            if variable == 0:
                subs_rate = r0
            else:
                subs_rate = r1
            astate = sa[variable]
            for bstate in range(4):
                sb = sa.copy()
                if astate == bstate:
                    continue
                sb[variable] = bstate
                idxb = mravel(sb)
                if (astate, bstate) in TRANSITIONS:
                    rate = subs_rate * kappa / (kappa + 2)
                else:
                    rate = subs_rate / (kappa + 2)
                Q[idxa, idxb] += rate

        # Add rates corresponding to synchronizations.
        # Synchronizations occur only between the second and third lineages.
        # Note that the rate is scaled by the the rate*time
        # of these lineages.
        for vsource, vsink in (1, 2), (2, 1):
            if sa[vsource] != sa[vsink]:
                sb = sa.copy()
                sb[vsink] = sa[vsource]
                idxb = mravel(sb)
                Q[idxa, idxb] += 0.5 * sync_rate * r1


def compute_prior_vector():
    """
    Compute the prior distribution over all 4*4*4 = 64 states.

    The prior distribution is uniform over all multivariate states
    for which all three variables share the same value in {0, 1, 2, 3}.
    """
    n = 4
    p = 1 / n
    prior = np.zeros(n*n*n, dtype=float)
    for state in product(range(n), repeat=3):
        if np.all(np.equal(state, state[0])):
            idx = mravel(state)
            prior[idx] = p
    return prior


def compute_rate_matrix(kappa, sync_rate, r0, r1):
    nstates = 4*4*4
    Q = np.zeros((nstates, nstates), dtype=float)
    add_rates(Q, kappa, sync_rate, r0, r1)
    exit_rates = Q.sum(axis=1)
    Q = Q - np.diag(exit_rates)
    return Q


def compute_likelihood_vector(prior_vector, rate_matrix):
    """
    Return a vector that gives the likelihood for each possible state.

    To compute the log likelihood given data and the parameterized rate matrix,
    compute the dot product of the observation count vector
    and the log of this likelihood vector.
    """
    P = expm(rate_matrix)
    likelihood_vector = prior_vector.dot(P)
    return likelihood_vector


def compute_p(prior_vector, X):
    """
    This is a helper function for computing the Fisher information matrix.

    """
    # Unpack the parameter values from the vector.
    kappa, sync_rate, r0, r1 = X
    rate_matrix = compute_rate_matrix(kappa, sync_rate, r0, r1)
    likelihood_vector = compute_likelihood_vector(prior_vector, rate_matrix)
    return likelihood_vector


def compute_kl_div(X, X_prime):
    """
    Compute K-L divergence.

    """
    prior_vector = compute_prior_vector()
    p_prime = compute_p(prior_vector, X_prime)
    p = compute_p(prior_vector, X)
    flat_div = kl_div(p_prime, p)
    return flat_div.sum()


def elementary(n, i):
    v = np.zeros(n, dtype=float)
    v[i] = 1
    return v


def compute_second_derivatives(f, delta, X):
    """
    Compute second and cross derivatives with finite differences.

    """
    n = X.shape[0]
    assert_equal(X.shape, (n,))
    M = np.zeros((n, n), dtype=float)
    delta2 = delta*delta
    for x in range(4):
        h = delta * elementary(n, x)
        M[x, x] = (f(X + h) - 2*f(X) + f(X - h)) / delta2
    for x in range(4):
        for y in range(x):
            h = delta * elementary(n, x)
            k = delta * elementary(n, y)
            value = (f(X+h+k) - f(X+h-k) - f(X-h+k) + f(X-h-k)) / (4*delta2)
            M[x, y] = value
            M[y, x] = value
    return M


def compute_fisher_information(X):
    """
    Return the Fisher information matrix using finite differences.

    The input X is a parameter vector as a length 4 ndarray.
    The output Fisher information matrix is a 4x4 matrix.
    If this matrix is positive definite, this would mean that the
    likelihood function is nicely behaved at the maximum likelihood estimate.

    """
    # Define the function whose partial derivatives are of interest.
    f = partial(compute_kl_div, X)

    # Estimate the second derivatives matrix.
    delta = 1e-4
    fisher_info = compute_second_derivatives(f, delta, X)
    return fisher_info


def run_analysis(kappa, sync_rate, r0, r1):
    """
    Look at the Fisher information at various points in the parameter space.

    For now use finite differences but maybe eventually
    use automatic differentiation.

    """
    X = np.array([kappa, sync_rate, r0, r1], dtype=float)
    print('kappa:', kappa)
    print('sync rate:', sync_rate)
    print('r0:', r0)
    print('r1:', r1)

    # Check the probability vector at the parameter values of interest.
    prior_vector = compute_prior_vector()
    p = compute_p(prior_vector, X)
    #print('probability vector:')
    #for i, value in enumerate(p):
        #print(i, value, sep='\t')
    #print()

    # Compute the fisher information.

    print('parameter vector:')
    print(X)

    print('fisher information matrix:')
    fisher_info = compute_fisher_information(X)
    print(fisher_info)

    print('spectrum of fisher information matrix:')
    print(eigvalsh(fisher_info))

    print('inverse of fisher information matrix:')
    print(inv(fisher_info))


def main():
    # The order is as follows:
    # kappa, sync_rate, r0, r1
    run_analysis(4, 2, 0.1, 0.2)
    run_analysis(4, 0.2, 0.1, 0.2)
    run_analysis(0.1, 1.0, 0.1, 0.2)
    run_analysis(4.0, 1.0, 0.001, 0.2)

if __name__ == '__main__':
    main()
