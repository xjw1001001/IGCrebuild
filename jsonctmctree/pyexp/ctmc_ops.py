"""
Define some custom linear operators.

Some of the operator norms of some operators can be computed efficiently.
This should really use more advanced abstract linear operator machinery.

"""
from __future__ import division, print_function, absolute_import

import numpy as np


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
