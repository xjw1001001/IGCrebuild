"""
A more advanced implementation of interface.py.

This implementation is more complicated and should be more efficient
than the naive implementation.

Below are some notes about which arrays are needed to meet which requests.
The set of resources is:
{
    low-memory conditional arrays,
    high-memory conditional arrays,
    high-memory subtree arrays,
}

The set of requirements is:

logl:
    ANY of the three resources
deri:
    high-memory conditional arrays
dwel:
    high memory subtree arrays
tran:
    high memory subtree arrays
root:
    ANY of the three resources
node:
    high-memory subtree arrays

"""
from __future__ import division, print_function, absolute_import

import sys

import numpy as np
from numpy.testing import assert_equal

from .expm_helpers import (
        ActionExpm,
        ImplicitDwellExpmFrechet,
        ImplicitTransitionExpmFrechetEx,
        )
from .common_likelihood import (
        get_conditional_likelihoods, get_subtree_likelihoods)
from .common_unpacking_ex import TopLevel, interpret_tree, interpret_root_prior
from .util import sparse_reduction
from . import expect
from . import ll
from .impl_naive import (
        _eagerly_precompute_dwell_objects,
        _apply_eagerly_precomputed_dwell_objects,
        _compute_transition_expectations,
        )


class InfeasibilityError(Exception):
    pass


class Reactor(object):
    """
    This is like a state machine.

    """
    def __init__(self, scene):
        self.scene = scene
        self.prior_distn = interpret_root_prior(scene)
        (
                self.T,
                self.root,
                self.edges,
                self.edge_rate_pairs,
                self.edge_process_pairs,
                ) = interpret_tree(scene)
        self.checked_feasibility = False
        self.likelihoods = None
        self.log_likelihoods = None
        self.node_to_subtree_likelihoods = None
        self.node_to_conditional_likelihoods = None
        self.node_to_marginal_distn = None

        # For each process, precompute the objects that are capable
        # of computing expm_mul and rate_mul for log likelihoods
        # and for its derivative with respect to edge-specific rates.
        self.expm_objects = []
        for p in scene.process_definitions:
            obj = ActionExpm(
                    scene.state_space_shape,
                    p.row_states,
                    p.column_states,
                    p.transition_rates)
            self.expm_objects.append(obj)


    def _delete_likelihoods(self, unmet_core_requests):
        if self.likelihoods is None:
            return False
        if not self.checked_feasibility:
            return False
        if unmet_core_requests & {'logl', 'root'}:
            return False
        self.likelihoods = None
        return True

    def _delete_log_likelihoods(self, unmet_core_requests):
        if self.log_likelihoods is None:
            return False
        if unmet_core_requests & {'logl'}:
            return False
        self.log_likelihoods = None
        return True

    def _delete_node_to_subtree_likelihoods(self, unmet_core_requests):
        if self.node_to_subtree_likelihoods is None:
            return False
        if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
            return False
        self.node_to_subtree_likelihoods = None
        return True

    def _delete_node_to_conditional_likelihoods(self, unmet_core_requests):
        if self.node_to_conditional_likelihoods is None:
            return False
        if unmet_core_requests & {'logl', 'deri', 'root'}:
            return False
        self.node_to_conditional_likelihoods = None
        return True

    def _delete_node_to_marginal_distn(self, unmet_core_requests):
        if self.node_to_marginal_distn is None:
            return False
        if unmet_core_requests & {'dwel', 'tran', 'node'}:
            return False
        self.node_to_marginal_distn = None
        return True


    def _check_feasibility(self):
        if self.checked_feasibility:
            return False
        if self.likelihoods is None:
            return False
        self.checked_feasibility = True
        if not np.all(self.likelihoods):
            raise InfeasibilityError
        return True


    def _create_log_likelihoods(self, unmet_core_requests):
        # Return True if an action was taken.
        if self.likelihoods is None:
            return False
        if 'logl' not in unmet_core_requests:
            return False
        self.log_likelihoods = np.log(self.likelihoods)
        return True

    #FIXME
    # site, edge, state
    #{D,S,W}NNLOGL : 3
    #{D,S,W}DNDERI : 3
    #{D,S,W}{D,W}{D,W}DWEL : 12
    #{D,S,W}{D,S,W}NTRAN : 9
    #{D,S,W}N{D,W}ROOT : 6
    #{D,S,W}N{D,W}NODE : 6

    def _respond_to_root(self, unmet_core_requests, requests, responses):
        if 'root' not in unmet_core_requests:
            return False
        if self.node_to_marginal_distn is not None:
        if self.likelihoods is None:
            return False
        self.self.root

        for i, request in enumerate(requests):
            if request.endswith('root'):

        return True

    def _respond_to_logl(self, unmet_core_requests, requests, responses):
        if 'logl' not in unmet_core_requests:
            return False
        return True

    def _respond_to_deri(self, unmet_core_requests, requests, responses):
        if 'deri' not in unmet_core_requests:
            return False
        return True

    def _respond_to_node(self, unmet_core_requests, requests, responses):
        if 'node' not in unmet_core_requests:
            return False
        return True

    def _respond_to_dwel(self, unmet_core_requests, requests, responses):
        if 'dwel' not in unmet_core_requests:
            return False
        return True

    def _respond_to_tran(self, unmet_core_requests, requests, responses):
        if 'tran' not in unmet_core_requests:
            return False
        return True


    def react(self, requests, responses):
        """
        This is called repeatedly, with some progress made in each call.
        The two input lists should have equal lengths, and this function
        should eventually replace each None entry of the responses list
        with a useful entry after enough calls.

        """
        # Get the set of unmet core properties.
        unmet_core_requests = set()
        for request, response in zip(requests, responses):
            if response is None:
                unmet_core_requests.add(request.propery[-4:])

        # Delete stuff.
        if self._gc_likelihoods(unmet_core_requests):
            return
        if self._gc_log_likelihoods(unmet_core_requests):
            return
        if self._gc_node_to_subtree_likelihoods(unmet_core_requests):
            return
        if self._gc_node_to_conditional_likelihoods(unmet_core_requests):
            return
        if self._gc_node_to_marginal_distn(unmet_core_requests):
            return

        # Check feasibility.
        if self._check_feasibility():
            return

        # Create intermediate arrays.

        # Respond to requests.
        if self._respond_to_root(unmet_core_requests, requests, responses):
            return
        if self._respond_to_logl(unmet_core_requests, requests, responses):
            return
        if self._respond_to_deri(unmet_core_requests, requests, responses):
            return
        if self._respond_to_node(unmet_core_requests, requests, responses):
            return
        if self._respond_to_dwel(unmet_core_requests, requests, responses):
            return
        if self._respond_to_tran(unmet_core_requests, requests, responses):
            return

        # Respond to log likelihoods requests.


        # Create the likelihoods array
        # from the node_to_subtree_likelihoods map.
        if self.node_to_subtree_likelihoods is not None:
            if self.likelihoods is None:
                if ('logl' in unmet_core_requests or
                        not self.checked_feasibility):
                    arr = self.node_to_subtree_likelihoods[self.root]
                    self.likelihoods = self.prior_distn.dot(arr)
                    assert_equal(len(self.likelihoods.shape), 1)
                    return

        # Create the likelihoods array
        # from the node_to_conditional_likelihoods map.
        if self.node_to_conditional_likelihoods is not None:
            if self.likelihoods is None:
                if ('logl' in unmet_core_requests or
                        not self.checked_feasibility):
                    arr = self.node_to_conditional_likelihoods[self.root]
                    self.likelihoods = self.prior_distn.dot(arr)
                    assert_equal(len(self.likelihoods.shape), 1)
                    return

        # Create the node_to_conditional_likelihoods map.
        if 

        if (not self.checked_feasibility and


            if self.likelihoods is None:
                if self.node_to_subtree_likelihoods is not None:
                    arr = self.node_to_subtree_likelihoods[self.root]
                    self.likelihoods = self.prior_distn.dot(arr)
                    assert_equal(len(self.likelihoods.shape), 1)
                    return
                if self.node_to_conditional_likelihoods is not None:
                    arr = self.node_to_conditional_likelihoods[self.root]
                    self.likelihoods = self.prior_distn.dot(arr)
                    return

                if not np.all(self.likelihoods):
                    return dict(
                            status = 'infeasible',
                            responses = None)

            if self.likelihoods is not None:


    # Compute an array related to expectations if requested.
    if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
        if debug:
            print('computing subtree likelihoods...', file=sys.stderr)
        store_all = True
        node_to_subtree_likelihoods = get_subtree_likelihoods(
                expm_objects,
                store_all,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                toplevel.scene.state_space_shape,
                toplevel.scene.observed_data.nodes,
                toplevel.scene.observed_data.variables,
                toplevel.scene.observed_data.iid_observations,
                )

    # Attempt to compute likelihoods
    # for the purpose of checking feasibility before progressing further.
    if likelihoods is None:
        if node_to_subtree_likelihoods is not None:
            arr = node_to_subtree_likelihoods[root]
            likelihoods = prior_distn.dot(arr)
            assert_equal(len(likelihoods.shape), 1)
            if not np.all(likelihoods):
                return dict(
                        status = 'infeasible',
                        responses = None)

    # Compute marginal distributions if necessary.
    if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
        if debug:
            print('computing marginal distributions...', file=sys.stderr)
        node_to_marginal_distn = expect.get_node_to_marginal_distn(
                expm_objects,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                toplevel.scene.state_space_shape,
                toplevel.scene.observed_data.nodes,
                toplevel.scene.observed_data.variables,
                toplevel.scene.observed_data.iid_observations,
                debug=debug)
        arr = []
        for i in range(toplevel.scene.node_count):
            arr.append(node_to_marginal_distn[i])
        full_node_array = np.array(arr).T

    def main(self, requests):
        responses = [None] * len(requests)
        try:
            while None in responses:
                self.react(requests, reponses)
            j_out = dict(
                    status = 'feasible',
                    responses = responses)
        except InfeasibilityError as e:
            j_out = dict(
                    status = 'infeasible',
                    responses = None)
        return j_out


def process_json_in(j_in, debug=False):
    """

    """
    toplevel = TopLevel(j_in)
    return Reactor(toplevel.scene).main(toplevel.requests)

    # Prepare a list of responses to be filled.
    # These are not necessarily filled in order from first to last.
    responses = [None] * len(requests)

    # Precompute the size of the state space.
    nstates = np.prod(scene.state_space_shape)
    iid_observation_count = scene.observed_data.iid_observations.shape[0]

    # Interpret the prior distribution by converting it to a dense array.
    prior_distn = interpret_root_prior(scene)

    # Interpret stuff related to the tree and its edges.
    info = interpret_tree(toplevel.scene)
    T, root, edges, edge_rate_pairs, edge_process_pairs = info
    nedges = len(edges)

    # For each process, precompute the objects that are capable
    # of computing expm_mul and rate_mul for log likelihoods
    # and for its derivative with respect to edge-specific rates.
    expm_objects = []
    for p in scene.process_definitions:
        obj = ActionExpm(
                scene.state_space_shape,
                p.row_states,
                p.column_states,
                p.transition_rates)
        expm_objects.append(obj)

    # Note that we haven't computed a few things yet.
    likelihoods = None
    log_likelihoods = None
    node_to_subtree_likelihoods = None
    node_to_marginal_distn = None
    full_node_array = None

    # Get the set of requested core properties.
    # Decisions will be made based on this set.
    unmet_core_requests = set(req.property[-4:] for req in requests)

    # Compute an array related to expectations if requested.
    if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
        if debug:
            print('computing subtree likelihoods...', file=sys.stderr)
        store_all = True
        node_to_subtree_likelihoods = get_subtree_likelihoods(
                expm_objects,
                store_all,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                toplevel.scene.state_space_shape,
                toplevel.scene.observed_data.nodes,
                toplevel.scene.observed_data.variables,
                toplevel.scene.observed_data.iid_observations,
                )

    # Attempt to compute likelihoods
    # for the purpose of checking feasibility before progressing further.
    if likelihoods is None:
        if node_to_subtree_likelihoods is not None:
            arr = node_to_subtree_likelihoods[root]
            likelihoods = prior_distn.dot(arr)
            assert_equal(len(likelihoods.shape), 1)
            if not np.all(likelihoods):
                return dict(
                        status = 'infeasible',
                        responses = None)

    # Compute marginal distributions if necessary.
    if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
        if debug:
            print('computing marginal distributions...', file=sys.stderr)
        node_to_marginal_distn = expect.get_node_to_marginal_distn(
                expm_objects,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                toplevel.scene.state_space_shape,
                toplevel.scene.observed_data.nodes,
                toplevel.scene.observed_data.variables,
                toplevel.scene.observed_data.iid_observations,
                debug=debug)
        arr = []
        for i in range(toplevel.scene.node_count):
            arr.append(node_to_marginal_distn[i])
        full_node_array = np.array(arr).T

    # In this naive implementation, always store all intermediate arrays.
    store_all = True
    node_to_conditional_likelihoods = get_conditional_likelihoods(
            expm_objects,
            store_all,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations,
            )

    # Get likelihoods at the root.
    # These are passed to the derivatives procedure,
    # to help compute the per-site derivatives of the log likelihoods
    # with respect to the log of the edge-specific rate scaling parameters.
    arr = node_to_conditional_likelihoods[root]
    likelihoods = prior_distn.dot(arr)
    assert_equal(len(likelihoods.shape), 1)

    # If the likelihood at any site is zero
    # then report infeasibility for the scene.
    if not np.all(likelihoods):
        return dict(
                status = 'infeasible',
                responses = None)

    # Apply the prior distribution and take logs of the likelihoods.
    log_likelihoods = np.log(likelihoods)

    # Compute the marginal distributions.
    # The full node array will have shape (nsites, nstatates, nnodes).
    if debug:
        print('computing marginal distributions...', file=sys.stderr)
    node_to_marginal_distn = expect.get_node_to_marginal_distn(
            expm_objects,
            node_to_subtree_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations,
            debug=debug)
    arr = []
    for i in range(toplevel.scene.node_count):
        arr.append(node_to_marginal_distn[i])
    full_node_array = np.array(arr).T

    # Compute the derivative of the likelihood
    # with respect to each edge-specific rate scaling parameter.
    # In this naive implementation compute all edge derivatives.
    requested_derivative_edge_indices = set(range(nedges))
    ei_to_derivatives = ll.get_edge_derivatives(
            expm_objects, requested_derivative_edge_indices,
            node_to_conditional_likelihoods, prior_distn,
            T, root, edges, edge_rate_pairs, edge_process_pairs,
            toplevel.scene.state_space_shape,
            toplevel.scene.observed_data.nodes,
            toplevel.scene.observed_data.variables,
            toplevel.scene.observed_data.iid_observations)

    # Fill an array with all unreduced derivatives.
    derivatives = np.empty((iid_observation_count, nedges))
    for ei, der in ei_to_derivatives.items():
        derivatives[:, ei] = der / likelihoods

    # Precompute dwell objects without regard to the requests.
    # This will obviously be changed for the less naive implementation.
    all_dwell_objects = _eagerly_precompute_dwell_objects(toplevel.scene)

    # Apply the precomputed dwell objects,
    # creating an array like (nsites, nedges, nstates) after the transposition.
    # A less naive implementation would compute this less eagerly.
    arr = []
    for dwell_state_index in range(nstates):
        dwell_objects = all_dwell_objects[dwell_state_index]
        arr.append(_apply_eagerly_precomputed_dwell_objects(
                toplevel.scene,
                expm_objects, dwell_objects, node_to_marginal_distn,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs))
    full_dwell_array = np.array(arr).T

    # Respond to each request, using a 'scene' common to all requests.
    responses = []
    for req in toplevel.requests:
        prefix, suffix = req.property[:3], req.property[-4:]
        observation_code, edge_code, state_code = prefix

        if suffix == 'logl':
            out = log_likelihoods
        elif suffix == 'deri':
            out = derivatives
        elif suffix == 'dwel':
            out = full_dwell_array
        elif suffix == 'tran':
            expm_transition_objects = []
            for p in toplevel.scene.process_definitions:
                obj = ImplicitTransitionExpmFrechetEx(
                        toplevel.scene.state_space_shape,
                        p.row_states,
                        p.column_states,
                        p.transition_rates,
                        req.transition_reduction.row_states,
                        req.transition_reduction.column_states,
                        req.transition_reduction.weights,
                        )
                expm_transition_objects.append(obj)
            arr = _compute_transition_expectations(
                toplevel.scene,
                expm_objects, expm_transition_objects, node_to_marginal_distn,
                node_to_subtree_likelihoods, prior_distn,
                T, root, edges, edge_rate_pairs, edge_process_pairs,
                debug=debug)
            out = np.array(arr).T
        elif suffix == 'root':
            out = full_node_array[:, :, root]
        elif suffix == 'node':
            out = full_node_array

        # Prepare to apply the reductions.
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
            indices = np.ravel_multi_index(
                    req.state_reduction.states.T,
                    toplevel.scene.state_space_shape)
            weights = req.state_reduction.weights
            out = sparse_reduction(out, indices, weights, reduction_axis)

        # Convert the ndarray to a list.
        # If the response is a zero ndim ndarray, then this will
        # convert to a single number.
        # Otherwise if the response is not an ndarray
        # (for example if it is a Python float),
        # then an AttributeError will be raised and the response
        # will be unchanged.
        try:
            out = out.tolist()
        except AttributeError as e:
            pass

        # Append the response.
        responses.append(out)

    # Return the response.
    return dict(status='feasible', responses=responses)
