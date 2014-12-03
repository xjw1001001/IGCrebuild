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


explanation of the input
------------------------

TODO


output
------

.. literalinclude:: out00.json
   :language: json
   :linenos:


explanation of the output
-------------------------

The output includes a status which can be "feasible", "infeasible", or "error",
together with an array of responses to the requests.
Because we issued only the single request for the logarithm of the probability
of ending in the original state, the response array has only a single entry,
namely the value of this logarithm computed as roughly -1
with some small numerical error.
