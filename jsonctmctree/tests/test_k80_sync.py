"""
This uses a model isomorphic to K80.

The transition/transversion ratio is assumed to be known.
Say that it is 4.0.
Each dot in a scatterplot represents a different combination
of 'substitution rate' and 'sync rate'.
The location of each dot is defined by the combination
of observed transition probability and observed transversion probability.

"""
from __future__ import division, print_function, absolute_import

from functools import partial
from itertools import product, permutations

import numpy as np
from scipy.linalg import expm
import scipy.stats

KAPPA = 4
TRANSITIONS = set([(0, 1), (1, 0), (2, 3), (3, 2)])


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
    #tips = sns.load_dataset('tips')
    #print(tips)
    #return

    n = 24
    #step_rates = 20 * np.random.rand(10)
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
    #pal = {0:'lightred', 0.5:'blue', 1:'seagreen', 2:'gray'}
    result = sns.lmplot(
            'transition_diff_prob',
            'transversion_diff_prob',
            df,
            hue='sync_rate',
            #palette=pal,
            #palette=pal,
            fit_reg=False)
    print('plot result:')
    print(result)
    print('showing the plot...')
    #sns.plt.show()
    sns.plt.savefig('K80-synchro-walks.svg')

    # Draw a plot.
    # This is from the seaborn tutorials.
    """
    g = sns.FacetGrid(df, hue='sync_rate', palette=pal, size=5)
    g.map(
            #plt.scatter,
            'state_difference_mean',
            'state_difference_variance',
            s=50,
            alpha=0.7,
            linewidth=0.5,
            edgecolor='white')
    g.add_legend()
    """


if __name__ == '__main__':
    main()
