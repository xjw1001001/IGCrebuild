"""
These are compatibility functions which should eventually be removed.

"""
from __future__ import division, print_function, absolute_import


import numpy as np
from scipy import sparse
from scipy.sparse import isspmatrix

__all__ = ['exact_1_norm', 'exact_inf_norm', 'trace', 'ident_like']


def exact_1_norm(A):
    if isspmatrix(A):
        return max(abs(A).sum(axis=0).flat)
    else:
        return np.linalg.norm(A, 1)

def exact_inf_norm(A):
    if isspmatrix(A):
        return max(abs(A).sum(axis=1).flat)
    else:
        return np.linalg.norm(A, np.inf)

def trace(A):
    if isspmatrix(A):
        return A.diagonal().sum()
    else:
        return np.trace(A)


def ident_like(A):
    if isspmatrix(A):
        return sparse.construct.eye(A.shape[0], A.shape[1],
                dtype=A.dtype, format=A.format)
    else:
        return np.eye(A.shape[0], A.shape[1], dtype=A.dtype)
