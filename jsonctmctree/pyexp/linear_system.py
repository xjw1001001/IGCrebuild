"""

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_

from .ctmc_ops import Propagator, ExplicitPropagator, RdOperator


class LinearSystem(object):
    """
    For now this just provides some help with linear operators.

    This code tries to be clever about using abstract linear operators vs.
    explicit sparse or dense matrix representations.

    One of the available operators is just the instantaneous transition
    rate operator.
    A second operator is available as a function of a requested time t.
    This corresponds to the action of the matrix whose (i, j) entry
    gives the probability of ending in state j after time t conditional
    on starting at state i.

    The style hint determines whether to use a dense vs. an abstract
    linear operator.  Note that regardless of the style hint,
    the input should still be a scipy sparse matrix.

    """
    def __init__(self, R, style='auto'):

        # Input validation.
        assert_(style in {'auto', 'dense', 'abstract'})
        if len(R.shape) != 2:
            raise ValueError
        if R.shape[1] != R.shape[0]:
            raise ValueError

        # Determine whether to use a dense matrix or not.
        n = R.shape[0]
        if style == 'auto':
            auto_threshold = 100
            if n < auto_threshold:
                use_dense_matrix = True
            else:
                use_dense_matrix = False
        elif style == 'dense':
            use_dense_matrix = True
        elif style == 'abstract':
            use_dense_matrix = False
        else:
            raise ValueError

        # R : sparse square matrix of non-negative off-diagonal rates
        self.shape = R.shape
        self.dtype = R.dtype

        # Compute exit rates.
        exit_rates = R.sum(axis=1).A.ravel()

        # Determine whether to use abstract or explicit linear operators.
        if use_dense_matrix:
            self._Q = R.A - np.diag(exit_rates)
            self._P = ExplicitPropagator(self._Q)
        else:
            d = -exit_rates
            mu = np.mean(d)
            op = RdOperator(R, d - mu)
            self._P = Propagator(op, mu)
            self._Q = RdOperator(R, d)

    @property
    def instantaneous_operator(self):
        return self._Q

    @property
    def propagator(self):
        # This could be converted to an operator using the following:
        # op = MatrixExponential(P, t)
        return self._P
