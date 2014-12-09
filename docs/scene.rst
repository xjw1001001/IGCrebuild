

.. _scene:

scene
-----

This object is one of the two attributes of the
jsonctmctree input,
the other being the array of :ref:`requests`.

The relatively simple scene below
has been taken from :ref:`example_1`
and is followed by detailed explanations of the
types and meanings of its components.

.. literalinclude:: tut00/in00.json
   :language: json
   :emphasize-lines: 2-26
   :linenos:


.. _shapes_and_sizes:

shapes and sizes
^^^^^^^^^^^^^^^^

This section defines some shapes and sizes of parts of the scenario.

    scene.node_count : integer
        The number of nodes in the branching timeline.

    scene.process_count : integer
        The number of distinct stochastic processes in the model.

    scene.state_space_shape : 1d array of integers
        The sizes of the state space space of each variable
        in the multivariate process.


.. _tree:

branching structure
^^^^^^^^^^^^^^^^^^^

This section defines the branching structure of the timeline,
including edge-specific rates and indicating which
processes act along which edges.

    scene.tree.row_nodes : 1d array of integers
        This array has one element for each edge of the tree,
        and the value of the element is the node index
        at the endpoint of the edge towards the root of the tree.
        If the directed graph of the rooted tree were represented
        by a sparse matrix, this array would represent row indices.

    scene.tree.column_nodes : 1d array of integers
        This is the complementary array of node indices.
        If the directed graph of the rooted tree were represented
        by a sparse matrix, this array would represent column indices.

    scene.tree.edge_rate_scaling_factors : 1d array of numbers
        Each edge is associated with a rate scaling factor.

    scene.tree.edge_processes : 1d array of integers
        Different edges are allowed to evolve according to different
        stochastic processes.
        This array gives the index of the stochastic process
        assigned to each edge.


.. _root_prior:

root prior
^^^^^^^^^^

This section defines the prior state distribution
at the root of the branching structure.

    scene.root_prior.states : 2d array of integers
        An array of multivariate states with nonzero
        probability at the root.

    scene.root_prior.probabilities : 1d array of numbers
        The distribution over multivariate states
        that have nonzero probability.


.. _process_definitions:

process definitions
^^^^^^^^^^^^^^^^^^^

This section defines the stochastic processes
which may act along one or more edges of the branching process.

    scene.process_definitions : 1d array of process objects
        Each object defines a stochastic process along one or more edges.

    scene.process_definitions[0].row_states : 2d array of integers
        Each entry of the array is a multivariate state.
        If the instantaneous transition rates were represented
        in matrix form, each entry would be the row index of a rate
        that is allowed to be nonzero.

    scene.process_definitions[0].column_states : 2d array of integers
        Each entry of the array is a multivariate state.
        If the instantaneous transition rates were represented
        in matrix form, each entry would be the column index of a rate
        that is allowed to be nonzero.

    scene.process_definitions[0].transition_rates : 1d array of numbers
        For each of the allowed transitions,
        this array contains the instantaneous rate of the transition.


.. _observed_data:

observed data
^^^^^^^^^^^^^

This section defines partial or complete observations
of the multivariate states at nodes of the branching structure,
by indicating which variables of the multivariate process
are observable at which nodes and providing multiple
such independent and identically distributed joint observations.

    scene.observed_data.nodes : 1d array of integers
        Indices of observable nodes.
        If multiple components of the multivariate process are observable
        at a node, then the node will be represented multiple times
        in this array.

    scene.observed_data.variables : 1d array of integers
        Indices of components of the multivariate process
        observable at the nodes indicated in the above array.

    scene.observed_data.iid_observations : 2d array of integers
        Observed component states of the multivariate process.
