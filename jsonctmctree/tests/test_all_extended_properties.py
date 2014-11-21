r"""
Superficially test all of the extended properties.

For all properties use the same scene.
The simplest nontrivial bivariate state space will be used.
This has two coevolving binary variables.
The resulting (2*2)x(2*2) rate matrix will be controlled
by three parameters:
    a: mutation rate from 0 to 1
    b: mutation rate from 1 to 0
    x: coalescence rate

In mathematical notation the state space is {0, 1} x {0, 1}.
In the notation used by the library,
this means that the shape of the state space is [2, 2],
the number of states in the multivariate state space is 4,
and the number of variables is 2.

The rate matrix defining the bivariate process is as follows:
     00       01       10   11
00 -a-a        a        a    0
01  a+x -a-b-x-x        0  b+x
10  b+x        0 -a-b-x-x  a+x
11    0        b        b -b-b

A second process will be defined by a setting x to zero.
This is the hadamard sum of two univariate processes,
because the two variables evolve independently.
     00       01       10   11
00 -a-a        a        a    0
01    a     -a-b        0    b
10    b        0     -a-b    a
11    0        b        b -b-b

A third process will have a different sparsity structure, like a gray code.
     00       01       10   11
00   -1        1        0    0
01    0       -1        0    1
10    1        0       -1    0
11    0        0        1   -1

The initial state distribution can be arbitrary and not symmetric, say
00 : 0.25
01 : 0.25
10 : 0.5
11 : 0

The tree could arbitrarily be like

     0
    / \
  (0) (1)
  /     \
 1       2
        / \
      (2) (3)
      /     \
     3       4

with a root, four edges, and five nodes.

We could say that edge 0 is controlled by the
process that has no dependence between variables,
while the other three edges have a dependence.
Say that edge 3 is controlled by the gray code process.

The data could arbitrarily consist of observations of both variables
at nodes 1 and 3, and only the second variable at nodes 2 and 4,
and no observation at the root.
The observed data themselves will be arbitrary.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_allclose, assert_equal

from jsonctmctree import interface

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


def _assert_list_shape(nested_lists, shape):
    arr = np.array(nested_lists)
    assert_equal(arr.shape, shape)


def _process_request(r):
    j_out = interface.process_json_in(dict(scene=_get_scene(), requests=[r]))
    assert_equal(set(j_out), {'status', 'responses'})
    assert_equal(j_out['status'], 'feasible')
    assert_equal(len(j_out['responses']), 1)
    out = j_out['responses'][0]
    prefix = r['property'][:3].lower()
    suffix = r['property'][-4:].lower()
    if suffix == 'node':
        assert_equal(len(np.array(out).shape), prefix.count('d') + 1)
    else:
        assert_equal(len(np.array(out).shape), prefix.count('d'))
    return out


"""
    'requests' : [
        {
            'property' : 'snnlogl'
        },
        {
            'property' : 'wwwdwel',
            'observation_reduction' : {
                'observation_indices' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'edge_reduction' : {
                'edges' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'state_reduction' : {
                'states' : [[0, 0], [0, 1], [0, 2]],
                'weights' : [1, 1, 1]}
        },
        {
            'property' : 'wsntran',
            'observation_reduction' : {
                'observation_indices' : [1, 1, 1],
                'weights' : [1, 1, 1]},
            'transition_reduction' : {
                'row_states' : [[0, 0], [0, 1], [0, 2]],
                'column_states' : [[0, 1], [0, 2], [0, 0]],
                'weights' : [1, 1, 1]}
        }
        ]
"""

# hard-code some request info
_sites = [0, 1, 2, 4, 3, 2]
_site_weights = [0.1, 0.1, 0.2, 0.3, 0.5, 0.8]

def test_logl():
    # {D,S,W}NNLOGL : 3
    observation_count = 5
    dnn = _process_request(dict(property='dnnlogl'))
    snn = _process_request(dict(property='snnlogl'))
    wnn = _process_request(dict(property='wnnlogl',
        observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights,
            ),
        ))
    _assert_list_shape(dnn, (observation_count, ))
    assert_allclose(np.sum(dnn, axis=0), snn)
    assert_allclose(interface._sparse_reduction(
        dnn, _sites, _site_weights, 0), wnn)

def test_deri():
    # {D,S,W}DNDERI : 3
    observation_count = 5
    edge_count = 4
    ddn = _process_request(dict(property='ddnderi'))
    sdn = _process_request(dict(property='sdnderi'))
    wdn = _process_request(dict(property='wdnderi',
        observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights),
        ))
    _assert_list_shape(ddn, (observation_count, edge_count))
    _assert_list_shape(sdn, (edge_count, ))
    _assert_list_shape(wdn, (edge_count, ))
    assert_allclose(np.sum(ddn, axis=0), sdn)
    assert_allclose(interface._sparse_reduction(
        ddn, _sites, _site_weights, 0), wdn)

def test_dwel():
    # {D,S,W}{D,W}{D,W}DWEL : 12
    observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights)
    edge_reduction = dict(
            edges=[0, 3, 2],
            weights=[0.4, 0.5, 2.0])
    state_reduction = dict(
            states=[[0, 0], [0, 1], [1, 0]],
            weights=[3, 3, 3])
    ddd = _process_request(dict(property='ddddwel'))
    sdd = _process_request(dict(property='sdddwel'))
    wdd = _process_request(dict(property='wdddwel',
        observation_reduction = observation_reduction))
    dwd = _process_request(dict(property='dwddwel',
        edge_reduction = edge_reduction))
    swd = _process_request(dict(property='swddwel',
        edge_reduction = edge_reduction))
    wwd = _process_request(dict(property='wwddwel',
        edge_reduction = edge_reduction,
        observation_reduction = observation_reduction))
    dww = _process_request(dict(property='dwwdwel',
        edge_reduction = edge_reduction,
        state_reduction = state_reduction))
    sww = _process_request(dict(property='swwdwel',
        edge_reduction = edge_reduction,
        state_reduction = state_reduction))
    www = _process_request(dict(property='wwwdwel',
        edge_reduction = edge_reduction,
        state_reduction = state_reduction,
        observation_reduction = observation_reduction))
    ddw = _process_request(dict(property='ddwdwel',
        state_reduction = state_reduction))
    sdw = _process_request(dict(property='sdwdwel',
        state_reduction = state_reduction))
    wdw = _process_request(dict(property='wdwdwel',
        state_reduction = state_reduction,
        observation_reduction = observation_reduction))

def test_tran():
    # {D,S,W}{D,S,W}NTRAN : 9
    observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights)
    edge_reduction = dict(
            edges=[0, 3, 2],
            weights=[0.4, 0.5, 2.0])
    transition_reduction = dict(
        row_states = [[0, 0], [0, 1], [1, 0]],
        column_states = [[1, 1], [1, 1], [0, 1]],
        weights = [1, 2, 3])
    ddn = _process_request(dict(property='ddntran',
        transition_reduction = transition_reduction))
    sdn = _process_request(dict(property='sdntran',
        transition_reduction = transition_reduction))
    wdn = _process_request(dict(property='wdntran',
        transition_reduction = transition_reduction,
        observation_reduction = observation_reduction))
    dsn = _process_request(dict(property='dsntran',
        transition_reduction = transition_reduction))
    ssn = _process_request(dict(property='ssntran',
        transition_reduction = transition_reduction))
    wsn = _process_request(dict(property='wsntran',
        transition_reduction = transition_reduction,
        observation_reduction = observation_reduction))
    dwn = _process_request(dict(property='dwntran',
        transition_reduction = transition_reduction,
        edge_reduction = edge_reduction))
    swn = _process_request(dict(property='swntran',
        transition_reduction = transition_reduction,
        edge_reduction = edge_reduction))
    wwn = _process_request(dict(property='wwntran',
        transition_reduction = transition_reduction,
        edge_reduction = edge_reduction,
        observation_reduction = observation_reduction))

def test_root():
    # {D,S,W}N{D,W}ROOT : 6
    observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights)
    state_reduction = dict(
            states=[[0, 0], [0, 1], [1, 0]],
            weights=[3, 3, 3])
    dnd = _process_request(dict(property='dndroot'))
    snd = _process_request(dict(property='sndroot'))
    wnd = _process_request(dict(property='wndroot',
        observation_reduction = observation_reduction))
    dnw = _process_request(dict(property='dnwroot',
        state_reduction = state_reduction))
    snw = _process_request(dict(property='snwroot',
        state_reduction = state_reduction))
    wnw = _process_request(dict(property='wnwroot',
        state_reduction = state_reduction,
        observation_reduction = observation_reduction))

def test_node():
    # {D,S,W}N{D,W}NODE : 6
    observation_reduction = dict(
            observation_indices=_sites,
            weights=_site_weights)
    state_reduction = dict(
            states=[[0, 0], [0, 1], [1, 0]],
            weights=[3, 3, 3])
    dnd = _process_request(dict(property='dndnode'))
    snd = _process_request(dict(property='sndnode'))
    wnd = _process_request(dict(property='wndnode',
        observation_reduction = observation_reduction))
    dnw = _process_request(dict(property='dnwnode',
        state_reduction = state_reduction))
    snw = _process_request(dict(property='snwnode',
        state_reduction = state_reduction))
    wnw = _process_request(dict(property='wnwnode',
        state_reduction = state_reduction,
        observation_reduction = observation_reduction))
