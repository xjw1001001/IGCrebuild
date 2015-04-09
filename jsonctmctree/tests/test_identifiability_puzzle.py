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
from numpy.testing import assert_allclose
import scipy.linalg

from jsonctmctree.sampling import assert_square_matrix
from jsonctmctree.ll import process_json_in
from jsonctmctree.sampling import(
        sample_distn,
        sample_time_reversible_rate_matrix,
        sample_time_nonreversible_rate_matrix,
        _brute_force_equilibrium,
        _geneconvify,
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


def test_identifiability():
    n = 4
    state_space_shape = (n, n)
    nstates = np.prod(state_space_shape)
    bivariate_states = list(product(range(n), repeat=2))
    bivariate_state_pairs = list(permutations(bivariate_states, 2))
    nrepeats = 10
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
