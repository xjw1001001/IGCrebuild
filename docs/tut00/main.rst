first tutorial
==============

In this tutorial we look at a simple continuous-time process
in which a binary variable begins in state 1 at
one end of an interval of length 1, 
and which stochastically evolves with instantaneous rate 1 towards
absorbing state 0.
This particular process is also morbidly known as a "pure-death" process
(as opposed to a "pure-birth" or "birth-death" process).

If we are interested in the probability to stay
in the initial state 1, we can see that this is simply the probability
of zero transitions within the interval during which 1 transition is expected.
According to the Poisson distribution,
this probability of no change is :math:`\frac{1}{e}`,
whose logarithm is -1.

This model is simple enough that using jsonctmctree to compute
this probability logarithm is somewhat awkward and unnecessary,
but because this is a jsonctmctree tutorial, we will do it anyway!

input
-----

.. literalinclude:: in00.json
   :language: json
   :linenos:

This input in the json meta-format formally defines
the scenario described in the introduction,
including a request for the logarithm of the probability
of no transition having occurred over the interval.

In the next few sections we will explain
the various parts of this unnecessarily complicated construction.


input shapes and sizes
^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 3-5
   :linenos:

This section defines some shapes and sizes of parts of the scenario.

    scene.node_count : integer
        Number of nodes in the branching timeline.

        In our example we have one node at each endpoint of the interval.

    scene.process_count : integer
        Number of distinct stochastic processes in the model.

        In our example we have only a single edge.
        Because each edge is associated with only one stochastic process,
        the process count in our example is 1.

    scene.state_space_shape : 1d array of integers
        The sizes of the state space space of each variable
        in the multivariate process.

        Our multivariate process has only a single variable
        so it is really a univariate process.
        This variable is binary so the size of its state space is 2.


input branching structure
^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 6-11
   :linenos:

This section defines the branching structure,
including edge-specific rates and indicating which
processes act along which edges.

    scene.tree : object
        The definition of the branching structure of the timeline.

        In our example we have only an interval, but we can think of this
        interval as the only edge in a not very interesting rooted tree
        whose root is node 0 and whose node at the opposite endpoint
        of the interval is 1.

    scene.tree.row_nodes : 1d array of integers
        This array has one element for each edge of the tree,
        and the value of the element is the node index
        at the endpoint of the edge towards the root of the tree.
        If the directed graph of the rooted tree were represented
        by a sparse matrix, this array would represent row indices.

        In our example we have one edge and two nodes,
        and so this array contains only the root node index which is 0.

    scene.tree.column_nodes : 1d array of integers
        This is the complementary array of node indices.
        If the directed graph of the rooted tree were represented
        by a sparse matrix, this array would represent column indices.

        In our example this array contains the only non-root node index.

    scene.tree.edge_rate_scaling_factors : 1d array of numbers
        Each edge is associated with a rate scaling factor.

        In our example we have only a single edge,
        and we set its scaling factor to 1.

    scene.tree.edge_processes : 1d array of integers
        Different edges are allowed to evolve according to different
        stochastic processes.
        This array gives the index of the stochastic process
        assigned to each edge.

        In our example we have only a single edge and only a single
        stochastic process.
        Because we use 0-based indexing, this array contains
        only a single entry which is 0.


input process definitions
^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 16-20
   :linenos:

This section defines the stochastic processes
which may act along one or more edges of the branching process.

    scene.process_definitions : 1d array of process objects
        Each object defines a stochastic process along one or more edges.

        In our example we have only a single edge so we need
        to define only a single stochastic process.
        Even if we had more edges, we could still get away with having
        only one process definition as long as all edges
        share the same stochastic process.

    scene.process_definitions[0].row_states : 2d array of integers
        Each entry of the array is a multivariate state.
        If the instantaneous transition rates were represented
        in matrix form, each entry would be the row index of a rate
        that is allowed to be nonzero.

        In our example only one transition between multivariate states
        is possible -- from 'multivariate' state [1] to state [0] -- and
        this array contains only the first of these two states.

    scene.process_definitions[0].column_states : 2d array of integers
        Each entry of the array is a multivariate state.
        If the instantaneous transition rates were represented
        in matrix form, each entry would be the column index of a rate
        that is allowed to be nonzero.

    scene.process_definitions[0].transition_rates : 1d array of numbers
        For each of the allowed transitions,
        this array contains the instantaneous rate of the transition.

        In our example we have only a single allowed transition
        and its instantaneous rate is 1.


input observed data
^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 21-25
   :linenos:

This section defines partial or complete observations
of the multivariate states at nodes of the branching structure.

    scene.observed_data : object
        Indicates which variables of the multivariate process
        are observable at which nodes, and provides multiple
        such independent and identically distributed joint observations.

    scene.observed_data.nodes : 1d array of integers
        Indices of observable nodes.
        If multiple components of the multivariate process are observable
        at a node, then the node will be represented multiple times
        in this array.

        In our example, we only care about observations at the
        node corresponding to the final endpoint of the interval.

    scene.observed_indices.variables : 1d array of integers
        Indices of components of the multivariate process
        observable at the nodes indicated in the above array.

        Our example uses a univariate process so it has only a single
        component -- component 0 -- and we only care about observations
        at one of the two endpoint nodes.

    scene.observed_indices.iid_observations : 2d array of integers
        Observed component states of the multivariate process.

        In our example we care only about a single observation,
        consisting of only a single variable observed at a single node.
        The observed value is 1 representing an observed absence
        of having transitioned to state 0 from state 1.


input requests
^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 27
   :linenos:

This section defines the requested properties of the scenario.

    requests : 1d array of request objects
        Array of requested reductions of requested properties
        of the 'scene'.
        
        In our example we only want a single logarithm of a probability,
        this array contains only a single request object.

    requests[0].property : string
        A code representing a reduction of a requested property.
        The first three letters represent possible reductions
        over observations, edges, and states respectively, with
        "D" representing no reduction,
        "S" representing unweighted summation,
        "W" representing weighted summation, and
        "N" meaning that reduction is inapplicable or disallowed.
        The suffix consisting of the last four letters
        belongs to one of the six codes
        {"LOGL", "DERI", "DWEL", "TRAN", "ROOT", "NODE"}.

        In our example, we request the property code "SNNLOGL" 
        corresponding to the log likelihood summed over observations.


output
------

.. literalinclude:: out00.json
   :language: json
   :linenos:

The output is also in the json meta-format.

    status : string
        Indicates feasibility or error status.

        In our example the status does not indicate an error.

    responses : 1d array of responses
        A response for each request.

        In our example we made only a single request
        for the logarithm of the probability of no transition
        to state 0, so this array has only a single entry.
        The logarithm of the probability was correctly computed
        as -1 with some small numerical error.
