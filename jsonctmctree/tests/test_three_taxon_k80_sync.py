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
from scipy.linalg import expm

KAPPA = 4
TRANSITIONS = set([(0, 1), (1, 0), (2, 3), (3, 2)])


def mravel(indices):
    return np.ravel_multi_index(indices, (4, 4, 4))


def add_rates(Q, kappa, sync_rate, r0, r1):
    multivariate_states = list(product(range(4), repeat=3))
    for sa in multivariate_states:
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
                if astate = bstate:
                    continue
                sb[variable] = bstate
                idxb = mravel(sb)
                if (astate, bstate) in TRANSITIONS:
                    rate = subs_rate * kappa / (kappa + 2)
                else:
                    rate = subs_rate / (kappa + 2)
                Q[idxa, idxb] += rate

        # Add rates corresponding to synchronizations.
        # Synchronizations occur only on the second and third lineages.

        #TODO unfinished below here...

        # Add rates corresponding to substitutions along other lineages.
        for bj in range(4):
            bi, bk = ai, ak
            idxb = mravel((bi, bj, bk))
            if ai == bi:
                continue
            elif (ai, bi) in TRANSITIONS:
                rate = r0 * kappa / (kappa + 2)
            else:
                rate = r0 / (kappa + 2)
            Q[idxa, idxb] = rate



def _fill_transition(expected_rate, n, Q, ai, aj, bi, bj):
    idxa = ai*n + aj
    idxb = bi*n + bj
    rate = expected_rate * KAPPA / (KAPPA + 2)
    Q[idxa, idxb] += rate


def _fill_transversion(expected_rate, n, Q, ai, aj, bi, bj):
    idxa = ai*n + aj
    idxb = bi*n + bj
    rate = expected_rate / (KAPPA + 2)
    Q[idxa, idxb] += rate


def _fill_sync(sync_rate, n, Q, ai, aj, bi, bj):
    idxa = ai*n + aj
    idxb = bi*n + bj
    Q[idxa, idxb] += 0.5 * sync_rate


def create_rate_matrix(expected_rate, sync_rate, n):
    Q = np.zeros((n*n, n*n), dtype=float)
    fill_transition = partial(_fill_transition, expected_rate, n, Q)
    fill_transversion = partial(_fill_transversion, expected_rate, n, Q)
    fill_sync = partial(_fill_sync, sync_rate, n, Q)
    for ai in range(n):
        for aj in range(n):

            # Add first variable transition and transversions.
            bi = ai
            for bj in 0, 1, 2, 3:
                if aj == bj:
                    continue
                elif (aj, bj) in TRANSITIONS:
                    fill_transition(ai, aj, bi, bj)
                else:
                    fill_transversion(ai, aj, bi, bj)

            # Add second variable transitions and transversions.
            bj = aj
            for bi in 0, 1, 2, 3:
                if ai == bi:
                    continue
                elif (ai, bi) in TRANSITIONS:
                    fill_transition(ai, aj, bi, bj)
                else:
                    fill_transversion(ai, aj, bi, bj)

            # Synchronize to the first variable.
            if ai != aj:
                bi, bj = ai, ai
                fill_sync(ai, aj, bi, bj)

            # Synchronize to the second variable.
            if ai != aj:
                bi, bj = aj, aj
                fill_sync(ai, aj, bi, bj)

    exit_rates = Q.sum(axis=1)
    Q = Q - np.diag(exit_rates)
    return Q


def compute_histogram(expected_rate, sync_rate, n):
    # This histogram is over (identities, transitions, transversions).
    Q = create_rate_matrix(expected_rate, sync_rate, n)
    P = expm(Q)
    initial_state = np.zeros(n*n, dtype=float)
    initial_state[0] = 1
    final_distribution = initial_state.dot(P)

    # Summarize the distribution over differences.
    histogram = np.zeros(3, dtype=float)
    for i in range(n):
        for j in range(n):
            idx = i*n + j
            p = final_distribution[idx]
            if i == j:
                d = 0
            elif (i, j) in TRANSITIONS:
                d = 1
            else:
                d = 2
            histogram[d] += p

    return histogram


def main():

    # Do a bunch of stuff from the seaborn tutorials.
    import pandas as pd
    import seaborn as sns
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    sns.set(style='white')
    np.random.seed(sum(map(ord, 'axis_grids')))

    n = 24
    subs_rates = np.linspace(0.5, 8, 41)

    transition_prob_col = []
    transversion_prob_col = []
    sync_rate_col = []
    subs_rate_col = []

    for sync_rate in 0, 0.5, 1, 2:
        npairs = len(subs_rates)
        for subs_rate in subs_rates:
            h = compute_histogram(subs_rate, sync_rate, n)
            transition_prob_col.append(h[1])
            transversion_prob_col.append(h[2])
        sync_rate_col.extend([sync_rate] * npairs)
        subs_rate_col.extend(subs_rates.tolist())


    # Create a data frame.
    # First make a Python dict.
    d = {
            'transition_diff_prob' : transition_prob_col,
            'transversion_diff_prob' : transversion_prob_col,
            'sync_rate' : sync_rate_col,
            'subs_rate' : subs_rate_col,
            }
    df = pd.DataFrame(data=d)

    # http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/
    # quantitative_linear_models.html#plotting-different-linear-relationships
    print('plotting a thing...')
    result = sns.lmplot(
            'transition_diff_prob',
            'transversion_diff_prob',
            df,
            hue='sync_rate',
            fit_reg=False)
    print('plot result:')
    print(result)
    print('showing the plot...')
    sns.plt.savefig('K80-synchro-walks.svg')


if __name__ == '__main__':
    main()
