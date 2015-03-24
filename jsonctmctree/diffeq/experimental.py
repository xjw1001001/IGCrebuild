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

_MMAX = 55

# This table helps to compute bounds.
# They seem to have been difficult to calculate, involving symbolic
# manipulation of equations, followed by numerical root finding.
_THETA = {
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


def mmax_to_pmax(mmax):
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


def _condition_3_13(A_1_norm, t, n0, mmax, ell):
    """
    A helper function for the _expm_multiply_* functions.

    Parameters
    ----------
    A_1_norm : float
        The precomputed 1-norm of A.
    n0 : int
        Number of columns in the _expm_multiply_* B matrix.
    mmax : int
        A value related to a bound.
    ell : int
        The number of columns used in the 1-norm approximation.
        This is usually taken to be small, maybe between 1 and 5.

    Returns
    -------
    value : bool
        Indicates whether or not the condition has been met.

    Notes
    -----
    This is condition (3.13) in Al-Mohy and Higham (2011).

    """

    # This is the rhs of equation (3.12).
    p_max = _compute_p_max(m_max)
    a = 2 * ell * p_max * (p_max + 3)

    # Evaluate the condition (3.13).
    b = _theta[m_max] / float(n0 * m_max)
    return A_1_norm * t <= a * b


class IterationStash(object):
    """
    Stash information for computing iteration counts.

    """
    def __init__(self, A):
        self._d = {}
        self._alpha = {}
        self._mmax = _MMAX
        self._pmax = mmax_to_pmax(self._mmax)
        self._A_1_norm = _exact_1_norm(A)

        # Initialize the S matrix.
        shape = (self._pmax - 1, self._mmax)
        self._S = np.zeros(shape)
        for p in range(2, self._pmax+1):
            for m in range(p*(p-1)-1, self._mmax+1):
                if m in _THETA:
                    self._S[p, m] = self.alpha(p) / _THETA[m]

    def d(self, p):
        # This calculation requires estimating
        # a root of the one-norm of a power of the A matrix.
        if p not in self._d:
            self._d[p] = foo
        return self._d[p]

    def alpha(self, p):
        if p not in self._alpha:
            self._alpha[p] = max(self.d(p), self.d(p+1))
        return self._alpha[p]

    def fragment_3_1(self, n0, t, ell):
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
            for m, theta in _theta.items():
                s = int(np.ceil(onenorm / theta))
                triples.append((m*s, m, s))
            ms, m, s = min(triples)
        else:
            m, s = self.cmstar(t)
        return m, s

    def condition_3_13(self, n0, t, ell):
        # This is the rhs of equation (3.12).
        a = 2 * ell * self._pmax * (self._pmax + 3)

        # Evaluate the condition (3.13).
        b = _THETA[self._mmax] / float(n0 * self._mmax)
        return self._A_1_norm * t <= a * b

    def cmstar(self, t):
        """
        This is a generalization of fragment 3.1.

        "and then for each t obtain C_mstar(t*A) as the smallest nonzero element
        in the matrix ceil(abs(t)*S)*diag(1, 2, ..., mmax), where mstar
        is the column index of the smallest element."

        a b c   u        au bv cw
        d e f *   v    = du ev fw
        g h i       w    gu hv iw

        [ a b c ]             au bv cw
        [ d e f ] * [u v w] = du ev fw
        [ g h i ]             gu hv iw

        """
        # Prepare a matrix that is like S except with inf replacing 0.
        M = np.where(self.S == 0, np.inf, self.S)

        # Transform entries of the matrix.
        abst = np.abs(t)
        row = np.arange(1, self._mmax + 1)
        M = np.ceil((abst * M) * row).astype(int)

        # Get the row and column of the smallest element.
        r, c = np.unravel_index(np.argmin(M), M.shape)
        mstar = c
        s = int(M[r, c])
        return mstar, s
