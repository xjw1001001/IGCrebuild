"""
A more advanced implementation of interface.py.

This implementation is more complicated and should be more efficient
than the naive implementation.

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


class InfeasibilityError(Exception):
    pass


#TODO the naive implementation should share this code
def apply_reductions(code, req, out):

    # Unpack the code.
    observation_code, edge_code, state_code = req

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
        indices = np.ravel_multi_index(
                req.state_reduction.states.T,
                toplevel.scene.state_space_shape)
        weights = req.state_reduction.weights
        out = sparse_reduction(out, indices, weights, reduction_axis)

    return out


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
        self.node_to_subtree_likelihoods = None
        self.node_to_conditional_likelihoods = None
        self.node_to_marginal_distn = None
        self.likelihoods = None
        self.log_likelihoods = None
        # If only the likelihoods and log likelihoods are required,
        # for example to check feasibility or to return log likelihoods
        # or to compute the posterior distribution at the root,
        # then only the root conditional likelihoods are required.
        self.root_conditional_likelihoods = None

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
        if unmet_core_requests & {'logl', 'deri'}:
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
        if unmet_core_requests & {'logl', 'deri'}:
            return False
        self.node_to_conditional_likelihoods = None
        return True

    def _delete_node_to_marginal_distn(self, unmet_core_requests):
        if self.node_to_marginal_distn is None:
            return False
        if unmet_core_requests & {'dwel', 'tran', 'root', 'node'}:
            return False
        self.node_to_marginal_distn = None
        return True

    def _delete_derivatives(self, unmet_core_requests):
        if self.derivatives is None:
            return False
        if unmet_core_requests & {'deri'}:
            return False
        self.derivatives = None
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


    def _create_likelihoods(self, unmet_core_requests):
        if self.likelihoods is not None:
            return False
        if self.checked_feasibility:
            if not (unmet_core_requests & {'logl', 'deri'}):
                return False
        if self.node_to_subtree_likelihoods is not None:
            arr = node_to_subtree_likelihoods[self.root]
        elif self.node_to_conditional_likelihoods is not None:
            arr = node_to_conditional_likelihoods[self.root]
        else:
            return False
        self.likelihoods = self.prior_distn.dot(arr)
        assert_equal(len(self.likelihoods.shape), 1)
        return True

    def _create_log_likelihoods(self, unmet_core_requests):
        if self.log_likelihoods is not None:
            return False
        if not (unmet_core_requests & {'logl'}):
            return False
        if self.likelihoods is None:
            return False
        self.log_likelihoods = np.log(self.likelihoods)
        return True

    def _create_derivatives(self, unmet_core_requests):
        if self.derivatives is not None:
            return False
        if not (unmet_core_requests & {'deri'}):
            return False
        if self.likelihoods is None:
            return False
        if self.node_to_conditional_likelihoods is None:
            return False

        # Compute the derivative of the likelihood
        # with respect to each edge-specific rate scaling parameter.
        #TODO do not necessarily request all edge derivatives
        nedges = len(self.edges)
        requested_derivative_edge_indices = set(range(nedges))
        ei_to_derivatives = ll.get_edge_derivatives(
                self.expm_objects, requested_derivative_edge_indices,
                self.node_to_conditional_likelihoods, self.prior_distn,
                self.T,
                self.root,
                self.edges,
                self.edge_rate_pairs,
                self.edge_process_pairs,
                self.scene.state_space_shape,
                self.scene.observed_data.nodes,
                self.scene.observed_data.variables,
                self.scene.observed_data.iid_observations)

        # Fill an array with all unreduced derivatives.
        iid_observation_count = len(self.scene.observed_data.iid_observations)
        self.derivatives = np.empty((iid_observation_count, nedges))
        for ei, der in ei_to_derivatives.items():
            self.derivatives[:, ei] = der / self.likelihoods
        return True

    def _create_node_to_conditional_likelihoods(self, unmet_core_requests):
        if self.node_to_conditional_likelihoods is not None:
            return False
        if self.checked_feasibility:
            if not (unmet_core_requests & {'logl', 'deri'}):
                return False
        #TODO restrict the requested number of arrays
        store_all = True
        self.node_to_conditional_likelihoods = get_conditional_likelihoods(
                self.expm_objects,
                store_all,
                self.T,
                self.root,
                self.edges,
                self.edge_rate_pairs,
                self.edge_process_pairs,
                self.scene.state_space_shape,
                self.scene.observed_data.nodes,
                self.scene.observed_data.variables,
                self.scene.observed_data.iid_observations,
                )
        return True

    def _create_node_to_subtree_likelihoods(self, unmet_core_requests):
        if self.node_to_subtree_likelihoods is not None:
            return False
        if self.checked_feasibility:
            if not (unmet_core_requests & {'dwel', 'tran', 'root', 'node'}):
                return False
        #TODO restrict the requested number of arrays
        store_all = True
        self.node_to_subtree_likelihoods = get_subtree_likelihoods(
                self.expm_objects,
                store_all,
                self.T,
                self.root,
                self.edges,
                self.edge_rate_pairs,
                self.edge_process_pairs,
                self.scene.state_space_shape,
                self.scene.observed_data.nodes,
                self.scene.observed_data.variables,
                self.scene.observed_data.iid_observations,
                )
        return True

    def _create_node_to_marginal_distn(self, unmet_core_requests):
        if self.node_to_marginal_distn is not None:
            return False
        if not (unmet_core_requests & {'dwel', 'tran', 'root', 'node'}):
            return False
        if node_to_subtree_likelihoods is None:
            return False
        debug = False
        self.node_to_marginal_distn = expect.get_node_to_marginal_distn(
                self.expm_objects,
                self.node_to_subtree_likelihoods, self.prior_distn,
                self.T,
                self.root,
                self.edges,
                self.edge_rate_pairs,
                self.edge_process_pairs,
                self.scene.state_space_shape,
                self.scene.observed_data.nodes,
                self.scene.observed_data.variables,
                self.scene.observed_data.iid_observations,
                debug=debug)
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
        if self.node_to_marginal_distn is None:
            return False
        full_array = self.node_to_marginal_distn[self.root].T
        for i, request in enumerate(requests):
            prefix = request.property[:3]
            suffix = request.property[-4:]
            if suffix == 'root':
                out = apply_reductions(prefix, request, full_array)
                responses[i] = out.tolist()
        return True

    def _respond_to_logl(self, unmet_core_requests, requests, responses):
        if 'logl' not in unmet_core_requests:
            return False
        if self.log_likelihoods is None:
            return False
        for i, request in enumerate(requests):
            prefix = request.property[:3]
            suffix = request.property[-4:]
            if suffix == 'logl':
                out = apply_reductions(prefix, request, self.log_likelihoods)
                responses[i] = out.tolist()
        return True

    def _respond_to_deri(self, unmet_core_requests, requests, responses):
        if 'deri' not in unmet_core_requests:
            return False
        if self.derivatives is None:
            return False
        for i, request in enumerate(requests):
            prefix = request.property[:3]
            suffix = request.property[-4:]
            if suffix == 'deri':
                out = apply_reductions(prefix, request, self.derivatives)
                responses[i] = out.tolist()
        return True

    def _respond_to_node(self, unmet_core_requests, requests, responses):
        if 'node' not in unmet_core_requests:
            return False
        if self.node_to_marginal_distn is None:
            return False
        arr = []
        for node in range(self.scene.node_count):
            arr.append(self.node_to_marginal_distn[node])
        full_node_array = np.array(arr).T
        for i, request in enumerate(requests):
            prefix = request.property[:3]
            suffix = request.property[-4:]
            if suffix == 'node':
                out = apply_reductions(prefix, request, full_node_array)
                responses[i] = out.tolist()
        return True

    #FIXME
    def _respond_to_dwel(self, unmet_core_requests, requests, responses):
        if 'dwel' not in unmet_core_requests:
            return False
        raise NotImplementedError
        return True

    #FIXME
    def _respond_to_tran(self, unmet_core_requests, requests, responses):
        if 'tran' not in unmet_core_requests:
            return False
        raise NotImplementedError
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
        if self._gc_node_to_conditional_likelihoods(unmet_core_requests):
            return
        if self._gc_node_to_subtree_likelihoods(unmet_core_requests):
            return
        if self._gc_node_to_marginal_distn(unmet_core_requests):
            return

        # Check feasibility.
        if self._check_feasibility():
            return

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

        # Create intermediate arrays.
        if self._create_likelihoods(unmet_core_requests):
            return
        if self._create_log_likelihoods(unmet_core_requests):
            return
        if self._create_conditional_likelihoods(unmet_core_requests):
            return
        if self._create_subtree_likelihoods(unmet_core_requests):
            return
        if self._create_derivatives(unmet_core_requests):
            return


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
    toplevel = TopLevel(j_in)
    return Reactor(toplevel.scene).main(toplevel.requests)
