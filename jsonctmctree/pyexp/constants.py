"""
Constants related to the matrix-exponential-vector product.

Al-Mohy and Higham
Computing the Action of the Matrix Exponential,
with an Application to Exponential Integrators (2011)
http://eprints.ma.man.ac.uk/1591/01/covered/MIMS_ep2010_30.pdf

Higham and Al-Mohy, 
Computing matrix functions (2010).
http://eprints.ma.man.ac.uk/1451/01/covered/MIMS_ep2010_18.pdf

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal

__all__ = ['MMAX', 'PMAX', 'THETA']


MMAX = 55

def _mmax_to_pmax(mmax):
    """
    Compute the largest positive integer p such that p*(p-1) <= m_max + 1.

    Parameters
    ----------
    mmax : int
        A count related to bounds.

    """
    sqrt_mmax = np.sqrt(mmax)
    p_low = int(np.floor(sqrt_mmax))
    p_high = int(np.ceil(sqrt_mmax + 1))
    pmax = max(p for p in range(p_low, p_high+1) if p*(p-1) <= mmax + 1)
    return pmax

PMAX = _mmax_to_pmax(MMAX)

assert_equal(PMAX, 8)


# This table helps to compute bounds.
# They seem to have been difficult to calculate, involving symbolic
# manipulation of equations, followed by numerical root finding.
THETA = {
        # The first 30 values are from table A.3 of Computing Matrix Functions.
        1: 2.29e-16,
        2: 2.58e-8,
        3: 1.39e-5,
        4: 3.40e-4,
        5: 2.40e-3,
        6: 9.07e-3,
        7: 2.38e-2,
        8: 5.00e-2,
        9: 8.96e-2,
        10: 1.44e-1,
        # 11
        11: 2.14e-1,
        12: 3.00e-1,
        13: 4.00e-1,
        14: 5.14e-1,
        15: 6.41e-1,
        16: 7.81e-1,
        17: 9.31e-1,
        18: 1.09,
        19: 1.26,
        20: 1.44,
        # 21
        21: 1.62,
        22: 1.82,
        23: 2.01,
        24: 2.22,
        25: 2.43,
        26: 2.64,
        27: 2.86,
        28: 3.08,
        29: 3.31,
        30: 3.54,
        # The rest are from table 3.1 of
        # Computing the Action of the Matrix Exponential.
        35: 4.7,
        40: 6.0,
        45: 7.2,
        50: 8.5,
        55: 9.9,
        }
