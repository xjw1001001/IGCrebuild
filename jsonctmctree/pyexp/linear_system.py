"""

"""
from __future__ import division, print_function, absolute_import

import numpy as np

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

    """
    def __init__(self, R):
        # R : sparse square matrix of non-negative off-diagonal rates
        self.shape = R.shape
        self.dtype = R.dtype

        # Compute exit rates.
        exit_rates = R.sum(axis=1).A.ravel()

        # Determine whether to use abstract or explicit linear operators.
        if exit_rates.shape[0] < 100:
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
