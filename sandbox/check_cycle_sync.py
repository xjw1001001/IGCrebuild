"""
"""
from __future__ import division, print_function, absolute_import

from functools import partial
from itertools import product, permutations

import numpy as np
from scipy.linalg import expm
import scipy.stats


def _fill_step(expected_rate, n, Q, ai, aj, bi, bj):
    idxa = ai*n + aj
    idxb = bi*n + bj
    Q[idxa, idxb] += 0.5 * expected_rate


def _fill_sync(sync_rate, n, Q, ai, aj, bi, bj):
    idxa = ai*n + aj
    idxb = bi*n + bj
    Q[idxa, idxb] += 0.5 * sync_rate


def create_rate_matrix(expected_rate, sync_rate, n):
    Q = np.zeros((n*n, n*n), dtype=float)
    fill_step = partial(_fill_step, expected_rate, n, Q)
    fill_sync = partial(_fill_sync, sync_rate, n, Q)
    for ai in range(n):
        for aj in range(n):

            # Add first variable transitions.
            bi = ai
            for bj in (aj-1)%n, (aj+1)%n:
                fill_step(ai, aj, bi, bj)

            # Add second variable transitions.
            bj = aj
            for bi in (ai-1)%n, (ai+1)%n:
                fill_step(ai, aj, bi, bj)

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
    Q = create_rate_matrix(expected_rate, sync_rate, n)
    P = expm(Q)
    initial_state = np.zeros(n*n, dtype=float)
    initial_state[0] = 1
    final_distribution = initial_state.dot(P)

    # Summarize the distribution over differences.
    distance_count = n // 2 + 1
    histogram = np.zeros(distance_count, dtype=float)
    for i in range(n):
        for j in range(n):
            idx = i*n + j
            p = final_distribution[idx]
            d = min((j-i)%n, (i-j)%n)
            histogram[d] += p

    return histogram


def get_moment_pairs(step_rates, sync_rate, n):
    distance_count = n // 2 + 1
    moment_pairs = []
    for step_rate in step_rates:
        histogram = compute_histogram(step_rate, sync_rate, n)
        xk = np.arange(distance_count)
        pk = histogram
        rv = scipy.stats.rv_discrete(values=(xk, pk))
        moment_pair = (rv.mean(), rv.var())
        moment_pairs.append(moment_pair)
    return moment_pairs


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
    step_rates = np.linspace(0.5, 8, 41)

    state_difference_mean_col = []
    state_difference_var_col = []
    sync_rate_col = []
    step_rate_col = []

    for sync_rate in 0, 0.5, 1, 2:
        moment_pairs = get_moment_pairs(step_rates, sync_rate, n)
        npairs = len(moment_pairs)
        means, variances = zip(*moment_pairs)

        state_difference_mean_col.extend(means)
        state_difference_var_col.extend(variances)
        sync_rate_col.extend([sync_rate] * npairs)
        step_rate_col.extend(step_rates.tolist())


    # Create a data frame.
    # First make a Python dict.
    d = {
            'state_difference_mean' : state_difference_mean_col,
            'state_difference_variance' : state_difference_var_col,
            'sync_rate' : sync_rate_col,
            'step_rate' : step_rate_col,
            }
    df = pd.DataFrame(data=d)

    # http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/
    # quantitative_linear_models.html#plotting-different-linear-relationships
    print('plotting a thing...')
    #pal = {0:'lightred', 0.5:'blue', 1:'seagreen', 2:'gray'}
    result = sns.lmplot(
            'state_difference_mean',
            'state_difference_variance',
            df,
            hue='sync_rate',
            #palette=pal,
            #palette=pal,
            fit_reg=False)
    print('plot result:')
    print(result)
    print('showing the plot...')
    #sns.plt.show()
    sns.plt.savefig('synchro-walks.svg')

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
