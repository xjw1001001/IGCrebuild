"""
These are some helper functions for HKY+IGC.

"""
from __future__ import print_function, division, absolute_import

import numpy as np
from scipy.special import logit, expit

__all__ = ['pack', 'unpack']


def _pack_acgt(pi):
    a, c, g, t = pi
    ag = a+g  # purines
    ct = c+t  # pyrimidines
    a_div_ag = a / ag
    c_div_ct = c / ct
    return logit([ag, a_div_ag, c_div_ct])


def _unpack_acgt(packed_acgt):
    ag, a_div_ag, c_div_ct = expit(packed_acgt)
    ct = 1 - ag
    a = a_div_ag * ag
    g = ag - a
    c = c_div_ct * ct
    t = ct - c
    return np.array([a, c, g, t])


def _pack_global_params(pi, kappa, tau):
    return np.concatenate([
        _pack_acgt(pi),
        np.log([kappa, tau])])


def _unpack_global_params(X):
    pi = _unpack_acgt(X[:3])
    kappa, tau = np.exp(X[3:])
    return pi, kappa, tau


def pack(distn, kappa, tau, rates):
    return np.concatenate((
        _pack_global_params(distn, kappa, tau),
        np.log(rates)))


def unpack(X):
    distn, kappa, tau = _unpack_global_params(X[:5])
    rates = np.exp(X[5:])
    return distn, kappa, tau, rates
