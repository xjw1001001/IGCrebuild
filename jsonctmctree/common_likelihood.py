"""
Code that is common for computing likelihood-related functions.

This is shared between the code that computes log likelihood and derivatives
and the code that computes posterior expectations of things like
linear combinations of labeled transitions on edges.

"""
from __future__ import division, print_function, absolute_import


def create_indicator_array(
        node,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations):
    """
    Create the initial array indicating observations.

    """
    nsites, nobservables = iid_observations.shape
    state_space_ndim = len(state_space_shape)
    state_space_axes = range(state_space_ndim)

    # Initialize the active array, initially in a high dimensional shape.
    # This array is large; for data with many iid sites,
    # such active arrays dominate the memory usage of the program.
    obs_shape = (nsites, ) + tuple(state_space_shape)
    obs = np.ones(obs_shape, dtype=float)

    # For each observable associated with the node under consideration,
    # apply the observation mask across all iid sites.
    local_observables = np.flatnonzero(observable_nodes == node)
    for idx in local_observables:
        states = iid_observations[:, idx]
        axis = observable_axes[idx]
        k = state_space_shape[axis]
        projection_shape = [k if i == axis else 1 for i in state_space_axes]
        mask_shape = (nsites, ) + tuple(projection_shape)
        ident = np.identity(k)
        obs *= np.take(ident, states, axis=0).reshape(mask_shape)

    # Reshape the observation array to 2d.
    # First collapse the dimensionality of the state space to 1d,
    # then transpose the array so that it has shape (nstates, nsites),
    # to prepare for P.dot(obs) where P has shape (nstates, nstates).
    obs = obs.reshape((nsites, np.prod(state_space_shape))).T

    # Return the observation indicator array.
    return obs


def get_conditional_likelihoods(
        # functions to compute expm_mul, per process (NEW)
        f,
        # boolean flag for storing all edges arrays (NEW)
        store_all,
        # original args
        T, root, edges, edge_rate_pairs, edge_process_pairs,
        state_space_shape,
        observable_nodes,
        observable_axes,
        iid_observations,
        ):
    """
    Recursively compute conditional likelihoods at the root.

    Attempt to order things intelligently to avoid using
    more memory than is necessary.

    The data provided by the caller gives us a sparse matrix
    of shape (nsites, nnodes, nstates).

    """
    child_to_edge = dict((tail, (head, tail)) for head, tail in edges)
    edge_to_rate = dict(edge_rate_pairs)
    edge_to_process = dict(edge_process_pairs)

    # For the few nodes that are active at a given point in the traversal,
    # we track a 2d array of shape (nsites, nstates).
    node_to_array = {}
    for node in get_node_evaluation_order(T, root):

        # When a node is activated, its associated array
        # is initialized to its observational likelihood array.
        arr = create_indicator_array(
                node,
                state_space_shape,
                observable_nodes,
                observable_axes,
                iid_observations)

        # When an internal node is activated,
        # this newly activated observational array is elementwise multiplied
        # by each of the active arrays of the child nodes.
        # The new elementwise product becomes the array
        # associated with the activated internal node.
        #
        # If we did not care about saving the per-node arrays,
        # then we could inactivate the child nodes and delete
        # their associated arrays, but because we want to re-use the
        # per-node arrays for edge length gradients, we keep them.
        for child in T.successors(node):
            arr *= node_to_array[child]
            if not store_all:
                del node_to_array[child]

        # When any node that is not the root is activated,
        # the matrix product P.dot(A) replaces A,
        # where A is the active array and P is the matrix exponential
        # associated with the parent edge.
        if node != root:
            edge = child_to_edge[node]
            edge_rate = edge_to_rate[edge]
            edge_process = edge_to_process[edge]
            arr = f[edge_process].expm_mul(edge_rate, arr)

        # Associate the array with the current node.
        node_to_array[node] = arr

    # If we had been deleting arrays as they become unnecessary for
    # the log likelihood calculation, then we would have only
    # a single active array remaining at this point, corresponding to the root.
    # But if we are saving the arrays for gradient calculations,
    # then we have more left.
    actual_keys = set(node_to_array)
    if store_all:
        desired_keys = set(T)
    else:
        desired_keys = {root}
    assert_equal(actual_keys, desired_keys)

    # Return the map from node to array.
    return node_to_array
