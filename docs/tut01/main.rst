example 2
=========

In this example we look at prior properties of the 2-state pure-death process.
In other words, we will look at properties of the model
when no data is observed.

.. note::

    Although we use some simple examples from survival theory
    or from the theory of stochastic birth-death processes,
    we do not want to suggest that jsonctmctree is only suited
    for these problems or even to suggest that it is especially
    well suited for this class of processes.
    It just so happened that its simplest possible applications
    fall into this class.


input
-----

.. literalinclude:: in00.json
   :language: json
   :linenos:

Except for the observed data section,
the 'scene' object used for this example is identical to that of example 1.

The next section will highlight the distinction between these two examples
in the way that the data is specified,
and the subsequent sections will highlight request-response pairs
for some prior summaries.


observed data
^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 21-25
   :linenos:

The highlighted section indicates that none of the variables
of the potentially multivariate process are observable at any of the nodes
of the potentially branching timeline,
and that we have a single 'observation' which contains
an empty array because nothing about the process is observable.


requests and responses
----------------------

In this section we highlight property requests
in the input json structure together with the corresponding
responses in the json output structure.


dwell proportions
^^^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 28-30
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :emphasize-lines: 4
   :linenos:

This is a request for "SDDDWEL" which is parsed as
the (S)ummation across observations,
for (D)istinct edges and (D)istinct states,
of the (DWEL)L proportions.
In this example, the value of ``responses[0][i][j]`` indicates
the proportion of time on edge i spent in state j.

The response values indicate that proportion
:math:`1 - \frac{1}{e}` of the time is expected to be spent in state 1.
This calculation agrees with the analytical calculation

.. math::

    & \int_0^1 \text{P(survival beyond time t)} \text{dt} \\
    &= \int_0^1 e^{-t} \text{dt} \\
    &= 1 - \frac{1}{e}


expected births
^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 31-38
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :emphasize-lines: 5
   :linenos:

This is a request for the "SSNTRAN" property which is parsed as
the (S)ummation across observations
of the (S)ummation across edges,
for which summation across states is (N)ot applicable,
of a linear combination of (TRAN)sition expectations.

The notation for a birth event is
a transition from multivariate state [0] to multivariate state [1].
Because we are requesting the expectation of only a single transition type,
the array of linear combination coefficients is simply [1].

As indicated by the name of the pure-death process,
we would not expect any births,
and the response of 0 agrees with this intuition.


expected deaths
^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 39-46
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :emphasize-lines: 6
   :linenos:

This similar "SSNTRAN" request asks for the expected number of death events.
Note that because of pecularities of the process,
this should be exactly equal to the probability
:math:`1 - \frac{1}{e}`
of survival over the interval, which is what we see.


states at nodes
^^^^^^^^^^^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 47-49
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :emphasize-lines: 7
   :linenos:

This is a request for the "SNDNODE" property which is parsed as
the (S)ummation across observations
for which summation across edges is (N)ot applicable, for (D)istinct states,
of the posterior probabilities at (NODE)s.
In this example, the value of ``responses[3][i][j]`` indicates
the posterior probability of state i at node j.

.. note::

    This (state, node) indexing is somewhat awkward.
    You might want each node to be associated with an array
    defining the posterior state distribution at that node,
    whereas in fact the response matrix has the transpose of the matrix.
    The reason for this awkwardness is to preserve consistency,
    so that arrays are indexed in the order (observation, edge, state)
    when possible.
    It would have been possible to add an fourth property prefix position
    indicating node reduction,
    but this addition would have its own interface complication tradeoffs.

We see that node 0
(also known as the initial endpoint of the interval, or the 'root' node)
has probability 0 of state 0
and probability 1 of state 1.
This is what we expect,
because we have specified this prior distribution
in the "scene.root_prior" section of the input json structure.

At node 1 representing the far endpoint of the interval,
we have probability :math:`1 - \frac{1}{e}` of state 0
and probability :math:`\frac{1}{e}` of state 1.
This is what we expect analytically.


conclusion
^^^^^^^^^^

Pretending for a moment that this example is a serious analysis
that misguidedly uses jsonctmctree instead of more competently using
the theory of survival analysis,
we could interpret the json responses as follows.

If we have a pure-death process initialized with state 1,
then after an elapsed time of 1 with instantaneous death rate of 1,
most of the time will be expected to have been spent in state 1,
although the most likely state after the elapsed time of 1 would
be state 0.
Furthermore, the expected number of death events is unsurprisingly
the same as the probability of not surviving beyond the end of the
time interval.
