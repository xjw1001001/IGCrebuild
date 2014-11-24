r"""
Compare outputs of different implementations of the interface.

See test_all_extended_properties.py for more information
about the scenario being tested.

"""
from __future__ import division, print_function, absolute_import

from itertools import product
import re

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from jsonctmctree import impl_naive, impl_v2
from jsonctmctree import common_unpacking_ex


#TODO share this across test files, in particular test_all_extended_properties
def _get_scene():
    a = 0.2
    b = 0.3
    x = 0.4
    return dict(
            node_count = 5,
            process_count = 2,
            state_space_shape = [2, 2],
            tree = dict(
                row_nodes = [0, 0, 2, 2],
                column_nodes = [1, 2, 3, 4],
                edge_rate_scaling_factors = [1.0, 2.0, 3.0, 4.0],
                edge_processes = [0, 1, 1, 2],
                ),
            root_prior = dict(
                states = [[0, 0], [0, 1], [1, 0]],
                probabilities = [0.25, 0.25, 0.5],
                ),
            process_definitions = [
                dict(
                    row_states = [
                        [0, 0], [0, 0], [0, 1], [0, 1],
                        [1, 0], [1, 0], [1, 1], [1, 1]],
                    column_states = [
                        [0, 1], [1, 0], [0, 0], [1, 1],
                        [0, 0], [1, 1], [0, 1], [1, 0]],
                    transition_rates = [a, a, a, b, b, a, b, b],
                    ),
                dict(
                    row_states = [
                        [0, 0], [0, 0], [0, 1], [0, 1],
                        [1, 0], [1, 0], [1, 1], [1, 1]],
                    column_states = [
                        [0, 1], [1, 0], [0, 0], [1, 1],
                        [0, 0], [1, 1], [0, 1], [1, 0]],
                    transition_rates = [a, a, a+x, b+x, b+x, a+x, b, b],
                    ),
                dict(
                    row_states = [[0, 0], [0, 1], [1, 1], [1, 0]],
                    column_states = [[0, 1], [1, 1], [1, 0], [0, 0]],
                    transition_rates = [1, 1, 1, 1],
                    ),
                ],
            observed_data = dict(
                nodes = [1, 1, 3, 3, 2, 4],
                variables = [0, 1, 0, 1, 1, 1],
                iid_observations = [
                    [0, 0, 0, 0, 0, 0],
                    [1, 1, 1, 1, 1, 1],
                    [0, 1, 0, 1, 0, 1],
                    [1, 1, 0, 0, 1, 1],
                    [1, 1, 1, 0, 0, 0],
                    ],
                ),
            )


def test_all_properties():

    # Initialize the scene dictionary.
    scene_dict = _get_scene()

    # hard-code all reductions
    observation_reduction = dict(
            observation_indices=[0, 1, 2, 4, 3, 2],
            weights=[0.1, 0.1, 0.2, 0.3, 0.5, 0.8])
    edge_reduction = dict(
            edges=[0, 3, 2],
            weights=[0.4, 0.5, 2.0])
    state_reduction = dict(
            states=[[0, 0], [0, 1], [1, 0]],
            weights=[3, 3, 3])
    transition_reduction = dict(
        row_states = [[0, 0], [0, 1], [1, 0]],
        column_states = [[1, 1], [1, 1], [0, 1]],
        weights = [1, 2, 3])

    # loop over all possible reductions
    #TODO implement all properties
    #core_properties = ('logl', 'deri', 'dwel', 'tran', 'root', 'node')
    core_properties = ('logl', 'deri', 'root', 'node')
    for components in product('dswn', 'dswn', 'dswn', core_properties):

        # Define the corresponding extended property.
        # Require that the extended property matches the regular expression.
        observation_code, edge_code, state_code, core_property = components
        extended_property = ''.join(components)
        if not re.match(common_unpacking_ex.request_regex, extended_property):
            continue

        print('requesting "%s"...' % extended_property)

        # Define the keyword arguments according to the extended property.
        # The 'tran' core property will always include
        # the 'transition_reduction' keyword argument.
        kwargs = dict()
        reductions = (observation_reduction, edge_reduction, state_reduction)
        names = ('observation_reduction', 'edge_reduction', 'state_reduction')
        codes = (observation_code, edge_code, state_code)
        for code, name, reduction in zip(codes, names, reductions):
            if code == 'w':
                kwargs[name] = reduction
        if core_property == 'tran':
            kwargs['transition_reduction'] = transition_reduction

        # Define the single request.
        request_dict = dict(property=extended_property)
        request_dict.update(kwargs)

        # Define the input structure including the scene and the single request.
        j_in = dict(
                scene=scene_dict,
                requests=[request_dict])

        # Compute the outputs using the two implementations.
        j_out_naive = impl_naive.process_json_in(j_in)
        j_out_v2 = impl_v2.process_json_in(j_in, debug=True)

        # Check some basic properties of the outputs.
        for j_out in j_out_naive, j_out_v2:
            assert_equal(set(j_out), {'status', 'responses'})
            assert_equal(j_out['status'], 'feasible')
            assert_equal(len(j_out['responses']), 1)

        # Compare the outputs for near equality.
        response_pairs = zip(j_out_naive['responses'], j_out_v2['responses'])
        for (response_naive, response_v2) in response_pairs:
            assert_allclose(response_naive, response_v2)

        print()
