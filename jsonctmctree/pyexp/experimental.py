"""
From the paper "Computing the Action of the Matrix Exponential":

If we wish to exponentiate the matrix t*A for several values of t then,
since alpha_p(t*A) = abs(t)*alpha_p(A), we can precompute the matrix S
with pmax - 1 rows and mmax columns given by
S_p_m = alpha_p(A) / theta_m when
2 <= p <= pmax and p*(p-1)-1 <= m <= mmax,
and 0 otherwise,
and then for each t obtain C_mstar(t*A) as the smallest nonzero element
in the matrix ceil(abs(t)*S)*diag(1, 2, ..., mmax), where mstar
is the column index of the smallest element.
Table 3.1 lists some of the theta_m values
corresponding to
u_s = tol = 2**-24 approx_eq 6e-8 (single precision) and
u_d = tol = 2**-53 approx_eq 1.1e-16 (double precision).
These values were determined as described in
'Computing Matrix Functions (2010)'

Higham and Al-Mohy, 
Computing matrix functions (2010).
http://eprints.ma.man.ac.uk/1451/01/covered/MIMS_ep2010_18.pdf

Al-Mohy and Higham
Computing the Action of the Matrix Exponential,
with an Application to Exponential Integrators (2011)
http://eprints.ma.man.ac.uk/1591/01/covered/MIMS_ep2010_30.pdf

"""
from __future__ import division, print_function, absolute_import

import numpy as np

from ._onenormest import onenormest
from .constants import MMAX, PMAX, THETA
from .sparse_dense_compat import exact_1_norm
from .basic_ops import PowerOperator


class RateMatrix(object):
    """
    This class is deliberately generally named.

    The input should be a scipy sparse matrix
    of non-negative off-diagonal rates.
    Diagonal entries should not be provided.

    """
    def __init__(self, A):
        self._A = A
        self._exit_rates = A.sum(axis=1).A.flatten()
        self._mean_exit_rate = np.mean(self._exit_rates)
        self._centered_exit_rates = self._exit_rates - self._mean_exit_rate

    def expm_multiply(self, t, B):
        pass

    def expm_adjoint_multiply(self, t, B):
        pass


class IterationStash(object):
    """
    Stash information for computing iteration counts.

    Some norms of matrix powers are pre-estimated for a given
    matrix with zero trace.

    These cached norms are subsequently used to compute
    iteration counts given a matrix scaling factor
    and the number of columns on the right hand side
    of the planned matrix-exponential-vector product.
    The first iteration count is related to the maximum
    degree of the matrix exponential Taylor expansion (although
    it may be stopped early depending on the data).
    The second iteration count is related to the number
    of requested segments that the interval should be broken into.

    """
    def __init__(self, A):
        self._A = A
        self._mmax = MMAX
        self._pmax = PMAX

        self._A_1_norm = A.one_norm()
        self._d = {1 : self._A_1_norm}
        self._alpha = {}

        # Initialize the S matrix.
        self._S = np.zeros((self._pmax+1, self._mmax+1))
        for p in range(2, self._pmax+1):
            for m in range(p*(p-1)-1, self._mmax+1):
                if m in THETA:
                    self._S[p, m] = self.alpha(p) / THETA[m]

        # Remove connections to some values that had been used
        # to create the _S matrix.
        # Some other values are kept.
        self._A = None
        self._d = None
        self._alpha = None

    def d(self, p):
        # This calculation requires computing a root of an estimate
        # of the one-norm of a power of the A matrix.
        if p not in self._d:
            op = PowerOperator(self._A, p)
            est = onenormest(op)
            self._d[p] = est ** (1 / p)
        return self._d[p]

    def alpha(self, p):
        if p not in self._alpha:
            self._alpha[p] = max(self.d(p), self.d(p+1))
        return self._alpha[p]

    def fragment_3_1(self, n0, t, ell=2):
        """
        Compute the order of the taylor expansion and the number of segments.

        If the matrix 1-norm is small then we only need the matrix 1-norm.
        Otherwise if the matrix 1-norm is large then we need to consider
        1-norm estimates of powers of the matrix.

        """
        if ell < 1:
            raise ValueError('expected ell to be a positive integer')
        elif self.condition_3_13(n0, t, ell):
            onenorm = self._A_1_norm * t
            triples = []
            for m, theta in THETA.items():
                s = int(np.ceil(onenorm / theta))
                triples.append((m*s, m, s))
            #print(triples)
            ms, m, s = min(triples)
        else:
            m, value = self.cmstar(t)
            s = max(int(value / m), 1)
        return m, s

    def condition_3_13(self, n0, t, ell):
        # This is the rhs of equation (3.12).
        a = 2 * ell * self._pmax * (self._pmax + 3)

        # Evaluate the condition (3.13).
        b = THETA[self._mmax] / float(n0 * self._mmax)
        return self._A_1_norm * t <= a * b

    def cmstar(self, t):
        """
        This is a generalization of fragment 3.1.

        "and then for each t obtain C_mstar(t*A) as the smallest nonzero element
        in the matrix ceil(abs(t)*S)*diag(1, 2, ..., mmax), where mstar
        is the column index of the smallest element."

        Traditional linear algebra notation
        -----------------------------------
        [ a b c ]   [ u     ]   [ au bv cw ]
        [ d e f ] * [   v   ] = [ du ev fw ]
        [ g h i ]   [     w ]   [ gu hv iw ]

        Numpy ndarray notation
        ----------------------
        [ a b c ]               [ au bv cw ]
        [ d e f ] * [ u v w ] = [ du ev fw ]
        [ g h i ]               [ gu hv iw ]

        """
        # Prepare a matrix that is like S except with inf replacing 0.
        M = np.where(self._S == 0, np.inf, self._S)

        # Transform entries of the matrix.
        abst = np.abs(t)
        row = np.arange(self._mmax + 1)
        row = np.where(row == 0, np.inf, row)
        M = np.ceil((abst * M) * row)

        # Get the row and value of the smallest element.
        r, c = np.unravel_index(np.argmin(M), M.shape)
        mstar = c
        value = int(M[r, c])
        return mstar, value
