"""
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_equal


__all__ = [
        'sparse_reduction',
        'apply_prefixed_reductions',
        'apply_reductions',
        ]


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


def apply_prefixed_reductions(state_space_shape, custom_prefix, req, out):
    """
    Apply reductions using a custom prefix.

    For example, if a reduction along one of the axes has already
    been performed using a more sophisticated method,
    then this function could be called with some letter other than
    {'d', 'n', 's', 'w'} at that position.

    """
    # Define the reduction codes.
    observation_code, edge_code, state_code = custom_prefix

    # Initialize the reduction axis.
    reduction_axis = 0

    # Apply the observation reduction if any.
    if observation_code == 'd':
        reduction_axis += 1
    elif observation_code == 's':
        out = np.sum(out, axis=reduction_axis)
    elif observation_code == 'w':
        indices = req.observation_reduction.observation_indices
        weights = req.observation_reduction.weights
        out = sparse_reduction(out, indices, weights, reduction_axis)

    # Apply the edge reduction if any.
    if edge_code == 'd':
        reduction_axis += 1
    elif edge_code == 's':
        out = np.sum(out, axis=reduction_axis)
    elif edge_code == 'w':
        indices = req.edge_reduction.edges
        weights = req.edge_reduction.weights
        out = sparse_reduction(out, indices, weights, reduction_axis)

    # Apply the state reduction if any.
    if state_code == 'd':
        reduction_axis += 1
    elif state_code == 's':
        out = np.sum(out, axis=reduction_axis)
    elif state_code == 'w':
        try:
            indices = np.ravel_multi_index(
                    req.state_reduction.states.T,
                    state_space_shape)
        except TypeError as e:
            raise Exception("%s; %s; %s" % (
                req.state_reduction.states.T,
                state_space_shape,
                e))
        weights = req.state_reduction.weights
        out = sparse_reduction(out, indices, weights, reduction_axis)

    return out


def apply_reductions(state_space_shape, req, out):
    """
    Apply reductions using the prefix defined by the request object.

    """
    # Unpack the prefix that defines the reduction axes.
    prefix = req.property[:3]
    return apply_prefixed_reductions(state_space_shape, prefix, req, out)
