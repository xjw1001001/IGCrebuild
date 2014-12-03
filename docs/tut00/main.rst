first tutorial
--------------

In this tutorial we look at a simple continuous-time process
in which a binary variable begins in state 1 at
one end of an interval of length 1, 
and which stochastically evolves with instantaneous rate 1 towards
absorbing state 0.

If we are interested in the probability to stay
in the initial state 1, we can see that this is simply the probability
of zero transitions within the interval during which 1 transition is expected.
According to the Poisson distribution,
this probability of no change is :math:`\frac{1}{e}`,
whose logarithm is :math:`-1`.

Because this is a jsonctmctree tutorial,
we will use it to compute this value!

.. literalinclude:: in00.json
   :language: json
   :linenos:

This logarithm is computed as roughly :math:`-1`
with some small numerical error.

.. literalinclude:: out00.json
   :language: json
   :linenos:
