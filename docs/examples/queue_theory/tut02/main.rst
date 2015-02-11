EM estimation
=============

In this example we use EM to estimate the two parameters
of a 3-state pure-death process.

Although EM is not built into jsonctmctree,
we can use it to compute expectations of the sufficient statistics,
and the maximization step is just a ratio.
We are continuing to use a pure-death model for our examples
because this process has the sparsest interesting rate matrix,
and is therefore in some sense the simplest interesting
continuous-time Markov process.


defining a model
----------------

In this section we will arbitrarily choose two
death rate parameter values

.. math::
   \mu_{2 \to 1} &= 0.25 \\
   \mu_{1 \to 0} &= 0.40

and use jsonctmctree to compute the distribution over states {0, 1, 2}
at the end of an elapsed time interval of length 1,
given that the initial state is 2.
We use methods that have been described more thoroughly in earlier examples.

.. note::

    The greek letter :math:`\mu` is used by convention in the birth-death
    process literature to indicate death rates.
    The corresponding birth rates are usually represented
    by :math:`\lambda` parameters.

input
^^^^^

.. literalinclude:: in00.json
   :language: json
   :emphasize-lines: 13, 19, 27
   :linenos:

The three highlighted lines of the input
show the initial state, the arbitrary transition rates,
and the requested property.
The requested property includes each
prior state probability at each node.

output
^^^^^^

.. literalinclude:: out00.json
   :language: json
   :linenos:

The output reassuringly shows that the initial (root) node
is unambiguously in state 2,
and it also shows the probability distribution
over states at the end of the unit time interval.
The interpretation of this distribution is that
after a time interval of length 1 has elapsed,
there is a 79% chance of remaining in state 2,
an 18% chance of having progressed to only state 1,
and a 4% chance of having progressed all the way to
the absorbing state 0.

The purpose of the remainder of this example is to show
how jsonctmctree can be used iteratively to
approximate the inverse of this calculation
(estimating the rate parameters from the final endpoint state distribution)
using a maximum likelihood inference strategy
called expectation maximization.
I should add the disclaimer that this is simply an example of
how jsonctmctree can be used to hack tiny models,
and it should not be taken as a recommendation to use this particular
solution for this particular problem.


estimating parameters
---------------------

In this section we pretend that we do not know
the arbitrary parameter values chosen above,
but we do know the distribution over
states after a unit interval of time has elapsed.


input scene
^^^^^^^^^^^

.. literalinclude:: in01.json
   :language: json
   :emphasize-lines: 19, 21-29
   :linenos:

To begin estimating parameters,
we construct a jsonctmctree input json structure
with the unknown transition parameters each arbitrarily initialized to 1
and with three observations at the endpoint node,
corresponding to the three possible endpoint states.


dwell expectations
^^^^^^^^^^^^^^^^^^

.. literalinclude:: in01.json
   :language: json
   :emphasize-lines: 32-38
   :linenos:

.. literalinclude:: out01.json
   :language: json
   :emphasize-lines: 4-8
   :linenos:

Our strategy to estimate each rate parameter will be to divide
the conditional expected number of transitions by the conditional expected
amount of opportunity to have made such a transition.
This section highlights the jsonctmctree request and response
related to the denominators of these estimates.

Notice that the request specifies
an :ref:`observation_reduction`
using weights corresponding to the distribution
over states at the final endpoint of the interval.

The response indicates that
about 88% of the elapsed time is expected to have been spent in state 2
and about 10% of the elapsed time is expected to have been spent in state 1.
In the next section, we will use these values as denominators
to produce the first EM iteration estimates.


transition expectations
^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in01.json
   :language: json
   :emphasize-lines: 39-61
   :linenos:

.. literalinclude:: out01.json
   :language: json
   :emphasize-lines: 9-10
   :linenos:

Here we show the same input/output pair of json structures,
but now the requests and responses related to conditional transition count
expectations have been highlighted.

The responses indicate that about 0.22 transitions from state 2 to state 1
are expected to have occurred conditional on the data,
and that about 0.04 transitions from state 1 to state 0
are expected to have occurred.

In fact for this particular model,
these conditional expectations do not depend on the parameter values,
given the probability distribution over states at the final endpoint.
So these requests are redundant.
We can see that the conditionally expected number of transitions from
state 2 to state 1 is equal to the probability of not ending up in
state 2, and the conditionally expected number of transitions from
state 1 to state 0 is equal to the probability of ending up in state 0.

Using these conditional transition count expectations together with the
conditional dwell time expectations, we can complete one
iteration of EM, updating our estimate
:math:`\hat{\mu}_{2 \to 1}` from its
arbitrary initial value of 1 to a more informed value of
0.221199216 / 0.8814781195026762 = 0.25094124.
Similarly we update
:math:`\hat{\mu}_{1 \to 0}` from its old value of 1
to a value of 0.040397988 / 0.10267733650267617 = 0.39344601
after this first iteration of expectation maximization.


second EM iteration
^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in02.json
   :language: json
   :emphasize-lines: 19
   :linenos:

.. literalinclude:: out02.json
   :language: json
   :linenos:

The second EM iteration uses almost the same input as the first iteration,
but using the updated rate parameters.

With these new conditional expectations we can update
:math:`\hat{\mu}_{2 \to 1}` from 0.25094124 to
0.22119921599999998 / 0.8846752438792551 = 0.2500343685780704
and update
:math:`\hat{\mu}_{1 \to 0}` from 0.39344601 to
0.040397988 / 0.10112290006727032 = 0.39949396203160625.

Note that the estimates 0.25003 and 0.39949
produced by this iteration seem to be approaching 0.25 and 0.40 respectively,
the original parameter values used to generate the
"observed" distribution over states at the final endpoint.

further EM iterations
^^^^^^^^^^^^^^^^^^^^^

This process could be repeated until some stopping criterion has been met.
According to Dempster, Laird, and Rubin (1977)
the sequence of parameter estimates :math:`x_k`
converges linearly to a limit :math:`L`
in the sense that

.. math::

   r = \lim_{k \to \infty}
   \dfrac{\mid x_{k+1} - L \mid}{\mid x_k - L \mid}

is a constant between 0 and 1 exclusive.

Of course if you wanted to do this for real,
you would write a program to do it instead of doing it by hand.
My hope is that this would not be such a formidable task,
because json readers and writers should be available
for most programming languages.
