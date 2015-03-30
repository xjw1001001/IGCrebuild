"""
Define some custom linear operators.

Some of the operator norms of some operators can be computed efficiently.
This should really use more advanced abstract linear operator machinery.

"""
from __future__ import division, print_function, absolute_import

import numpy as np

from .experimental import IterationStash


__all__ = ['RdOperator', 'RdcOperator', 'RdCOperator']


class _HighLevelInterface(object):
    # This interface expects to be mixed into a class with the following members
    # _matmat (function)
    # _transpose (function)
    # _adjoint (function)

    @property
    def T(self):
        return self._transpose()

    @property
    def H(self):
        return self._adjoint()

    def dot(self, other):
        return self._matmat(other)



class _ConcreteInterface(object):
    # This interface expects to be mixed into a class with the following members
    # abs_sum_axis_0 (function)
    # abs_sum_axis_1 (function)
    # _matmat (function)
    # _my_adjoint_matmat (function)
    # shape (property)
    # dtype (property)

    def _init_concrete_cache(self):
        self._one_norm = None
        self._inf_norm = None
        self.__adj = None

    def one_norm(self):
        if self._one_norm is None:
            self._one_norm = np.max(self.abs_sum_axis_0())
        return self._one_norm

    def inf_norm(self):
        if self._inf_norm is None:
            self._inf_norm = np.max(self.abs_sum_axis_1())
        return self._inf_norm

    def _transpose(self):
        if self.__adj is None:
            self.__adj = _ExtendedAdjointOperator(self)
        return self.__adj

    def _adjoint(self):
        return self._transpose()


class _VanillaAdjointOperator(_HighLevelInterface):
    # The forward operator is expected to meet a few extra conditions.
    # It should be real-valued (so the adjoint is the transpose).
    # It should also be square but not necessarily symmetric.
    # This class is not mixed with _ConcreteInterface.
    def __init__(self, L):
        self._L = L
        self.args = (L, )
    @property
    def shape(self):
        return self._L.shape
    @property
    def dtype(self):
        return self._L.dtype
    def _matmat(self, other):
        return self._L._my_adjoint_matmat(other)
    def _adjoint(self):
        return self._L
    def _transpose(self):
        return self._L


class _ExtendedAdjointOperator(_HighLevelInterface):
    # The forward operator is expected to meet a few extra conditions.
    # It should be real-valued (so the adjoint is the transpose),
    # and it should have easily computable row and column vector 1-norms.
    # It should also be square but not necessarily symmetric.
    # This class is not mixed with _ConcreteInterface.
    def __init__(self, L):
        self._L = L
        self.args = (L, )
    @property
    def shape(self):
        return self._L.shape
    @property
    def dtype(self):
        return self._L.dtype
    def _matmat(self, other):
        return self._L._my_adjoint_matmat(other)
    def abs_sum_axis_0(self):
        return self._L.abs_sum_axis_1()
    def abs_sum_axis_1(self):
        return self._L.abs_sum_axis_0()
    def one_norm(self):
        return self._L.inf_norm()
    def inf_norm(self):
        return self._L.one_norm()
    def _adjoint(self):
        return self._L
    def _transpose(self):
        return self._L


class PowerOperator(_HighLevelInterface):
    def __init__(self, L, p):
        self._L = L
        self._p = p
        self.__adj = None
    def _matmat(self, other):
        for i in range(self._p):
            other = self._L.dot(other)
        return other
    def _my_adjoint_matmat(self, other):
        for i in range(self._p):
            other = self._L.H.dot(other)
        return other
    def _transpose(self):
        if self.__adj is None:
            self.__adj = _VanillaAdjointOperator(self)
        return self.__adj
    def _adjoint(self):
        return self._transpose()


class _ExtendedMatrixOperator(_HighLevelInterface, _ConcreteInterface):
    # This wraps a sparse matrix and has an interface common to this module.
    # The matrix should be square and real-valued.
    def __init__(self, M):
        self.dtype = M.dtype
        self.shape = M.shape
        self._M = M
        self.args = (M, )
        self._MT = None
        self._M_abs = None
        self._abs_sum_axis_0 = None
        self._abs_sum_axis_1 = None
        self._init_concrete_cache()

    def abs_sum_axis_0(self):
        if self._abs_sum_axis_0 is None:
            if self._M_abs is None:
                self._M_abs = np.absolute(self._M)
            self._abs_sum_axis_0 = self._M_abs.sum(axis=0).A.ravel()
        return self._abs_sum_axis_0

    def abs_sum_axis_1(self):
        if self._abs_sum_axis_1 is None:
            if self._M_abs is None:
                self._M_abs = np.absolute(self._M)
            self._abs_sum_axis_1 = self._M_abs.sum(axis=1).A.ravel()
        return self._abs_sum_axis_1

    def _matmat(self, other):
        return self._M.dot(other)

    def _my_adjoint_matmat(self, other):
        if self._MT is None:
            self._MT = self._M.T
        return self._MT.dot(other)


class RdOperator(_HighLevelInterface, _ConcreteInterface):
    """
    This is a custom linear operator.

    It is the sum of a square sparse matrix R and a diagonal matrix d.
    The R component is assumed to have entrywise non-negative
    floating point values with zero diagonal.
    The d component is represented by a 1d ndarray
    with floating point entries.
    No complex numbers are involved,
    so conjugate transposes are just transposes.
    The 1-norm of this operator can be computed without much difficulty.

    """
    def __init__(self, R, d):
        self.dtype = R.dtype
        self.shape = R.shape
        self._R = R
        self._d = d
        self._d_abs = np.abs(d)
        self.args = R, d
        self._RT = None
        self._abs_sum_axis_0 = None
        self._abs_sum_axis_1 = None
        self._init_concrete_cache()

    def abs_sum_axis_0(self):
        if self._abs_sum_axis_0 is None:
            self._abs_sum_axis_0 = self._R.sum(axis=0).A.ravel() + self._d_abs
        return self._abs_sum_axis_0

    def abs_sum_axis_1(self):
        if self._abs_sum_axis_1 is None:
            self._abs_sum_axis_1 = self._R.sum(axis=1).A.ravel() + self._d_abs
        return self._abs_sum_axis_1

    def _matmat(self, other):
        return self._R.dot(other) + self._d[:, np.newaxis] * other

    def _my_adjoint_matmat(self, other):
        if self._RT is None:
            self._RT = self._R.T
        return self._RT.dot(other) + self._d[:, np.newaxis] * other


class RdcOperator(_HighLevelInterface, _ConcreteInterface):
    # R+d  c
    #  0  R+d
    def __init__(self, Rd, c):
        self.dtype = Rd.dtype
        self.shape = Rd.shape[0]*2, Rd.shape[1]*2
        self._Rd = Rd
        self._c = c
        self._c_abs = np.abs(c)
        self.args = Rd, c
        self._abs_sum_axis_0 = None
        self._abs_sum_axis_1 = None
        self._init_concrete_cache()

    def abs_sum_axis_0(self):
        if self._abs_sum_axis_0 is None:
            self._abs_sum_axis_0 = np.concatenate((
                self._Rd.abs_sum_axis_0(),
                self._Rd.abs_sum_axis_0() + self._c_abs))
        return self._abs_sum_axis_0

    def abs_sum_axis_1(self):
        if self._abs_sum_axis_1 is None:
            self._abs_sum_axis_1 = np.concatenate((
                self._Rd.abs_sum_axis_1() + self._c_abs,
                self._Rd.abs_sum_axis_1()))
        return self._abs_sum_axis_1

    def _matmat(self, other):
        n = self._c.size
        M = np.empty_like(other)
        M[:n, :] = (self._Rd.dot(other[:n, :]) +
                    self._c[:, np.newaxis] * other[n:, :])
        M[n:, :] = self._Rd.dot(other[n:, :])
        return M

    def _my_adjoint_matmat(self, other):
        n = self._c.size
        M = np.empty_like(other)
        M[:n, :] = self._Rd._my_adjoint_matmat(other[:n, :])
        M[n:, :] = (self._Rd._my_adjoint_matmat(other[n:, :]) +
                    self._c[:, np.newaxis] * other[:n, :])
        return M


class RdCOperator(_HighLevelInterface, _ConcreteInterface):
    # R+d  C
    #  0  R+d
    def __init__(self, Rd, C):
        self.dtype = Rd.dtype
        self.shape = Rd.shape[0]*2, Rd.shape[1]*2
        self._Rd = Rd
        self._C = _ExtendedMatrixOperator(C)
        self._CH = None
        self.args = Rd, C
        self._abs_sum_axis_0 = None
        self._abs_sum_axis_1 = None
        self._init_concrete_cache()

    def abs_sum_axis_0(self):
        if self._abs_sum_axis_0 is None:
            self._abs_sum_axis_0 = np.concatenate((
                self._Rd.abs_sum_axis_0(),
                self._Rd.abs_sum_axis_0() + self._C.abs_sum_axis_0()))
        return self._abs_sum_axis_0

    def abs_sum_axis_1(self):
        if self._abs_sum_axis_1 is None:
            self._abs_sum_axis_1 = np.concatenate((
                self._Rd.abs_sum_axis_1() + self._C.abs_sum_axis_1(),
                self._Rd.abs_sum_axis_1()))
        return self._abs_sum_axis_1

    def _matmat(self, other):
        n = self._C.shape[0]
        M = np.empty_like(other)
        M[:n, :] = (self._Rd.dot(other[:n, :]) +
                    self._C.dot(other[n:, :]))
        M[n:, :] = self._Rd.dot(other[n:, :])
        return M

    def _my_adjoint_matmat(self, other):
        n = self._C.shape[0]
        if self._CH is None:
            self._CH = self._C.H
        M = np.empty_like(other)
        M[:n, :] = self._Rd._my_adjoint_matmat(other[:n, :])
        M[n:, :] = (self._Rd._my_adjoint_matmat(other[n:, :]) +
                    self._CH.dot(other[:n, :]))
        return M


def _expm_product_helper(A, mu, iteration_stash, t, B):
    # Estimate expm(t*M).dot(B).
    # A = M - mu*I
    # mu = mean(trace(M))
    # The iteration stash helps to compute numbers of iterations to use.
    # t is a scaling factor.
    # B is the input matrix for the linear operator.

    # Compute some input-dependent constants.
    tol = np.ldexp(1, -53)
    n0 = B.shape[1]
    m, s = forward_iteration_stash.fragment_3_1(n0, t)

    # Get the lapack function for computing matrix norms.
    lange, = get_lapack_funcs(('lange',), (B,))

    F = B
    eta = np.exp(t*mu / float(s))
    for i in range(s):
        c1 = lange('i', B)
        for j in range(m):
            coeff = t / float(s*(j+1))
            B = coeff * A.dot(B)
            c2 = lange('i', B)
            F = F + B
            if c1 + c2 <= tol * lange('i', F):
                break
            c1 = c2
        F = eta * F
        B = F
    return F


def _Propagator(object):
    """
    Wraps a linear operator.

    Given a linear operator L and a scaling factor t,
    this provides the operator (and adjoint operator)
    corresponding to the matrix exponential of the scaled operator L*t.
    It is intended to be wrapped by a function
    that provides the corresponding operator specific to a given t.

    The point of this class is basically to cache properties
    that are the same across different values of t.

    """
    def __init__(self, A, mu):
        # A = M - mu*I is an abstract linear operator
        # whose 1-norm is directly accessible.
        # mu is the mean trace of M.
        self._A = A
        self._mu = mu
        self._forward_iteration_stash = None
        self._adjoint_iteration_stash = None

    def _parameterized_matvec(self, t, B):
        # Approximate expm(M*t).dot(B).
        # t is a scaling factor of L
        # B the input matrix of the linear function
        if self._forward_iteration_stash is None:
            self._forward_iteration_stash = IterationStash(self._A)
        return _expm_product_helper(
                self._A, self._mu, self._forward_iteration_stash, t, B)

    def _parameterized_adjoint_matvec(self, t, B):
        # Approximate expm(M.H*t).dot(B).
        # t is a scaling factor of L
        # B the input matrix of the adjoint linear function
        if self._adjoint_iteration_stash is None:
            self._adjoint_iteration_stash = IterationStash(self._A.H)
        return _expm_product_helper(
                self._A.H, self._mu, self._adjoint_iteration_stash, t, B)


def MatrixExponential(_HighLevelInterface):
    # The input is already a propagator; this just scales by a specific t.
    def __init__(self, P, t):
        # P is a propagator
        # t is the scaling factor
        self._P = P
        self._t = t
        self.__adj = None

    def _matvec(self, B):
        # B is the input matrix of the linear function.
        return self._P._parameterized_matvec(self._t, B)

    def _my_adjoint_matvec(self, B):
        # B is the input matrix of the adjoint linear function.
        return self._P._parameterized_adjoint_matvec(self._t, B)

    def _transpose(self):
        # This uses a vanilla adjoint because norms are unavailable.
        if self.__adj is None:
            self.__adj = _VanillaAdjointOperator(self)
        return self.__adj

    def _adjoint(self):
        # This is the same as the transpose because the entries are real.
        return self._transpose()
