"""
Utility functions.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal

__all__ = [
        'assert_square',
        'sparse_reduction']


def assert_square(M):
    assert_equal(len(Q.shape), 2)
    assert_equal(Q.shape[0], Q.shape[1])


def sparse_reduction(A, indices, weights, axis):
    """
    A sparse linear combination along an axis.

    Parameters
    ----------
    A : ndarray
        the ndarray to reduce
    indices : 1d ndarray
        the indices to be used in the reduction along the axis of interest
    weights : 1d ndarray
        the array of weights to be used in the reduction
    axis : int
        the axis along which to reduce the input ndarray

    """
    A = np.asarray(A)
    indices = np.asarray(indices)
    weights = np.asarray(weights)
    assert_equal(indices.ndim, 1)
    assert_equal(weights.ndim, 1)
    assert_equal(indices.shape, weights.shape)
    return np.tensordot(
            np.take(A, indices, axis=axis),
            weights,
            axes=([axis], [0]))
