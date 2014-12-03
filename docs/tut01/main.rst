first tutorial
--------------

In this tutorial we look at a univariate binary continuous-time
Markov process defined by instantaneous rates

.. math::

   \lambda_{0 \to 1} &= 1 \\
   \lambda_{1 \to 0} &= 1

along an interval of length 1.
The process is assumed to be at stationarity,
so the prior distribution over states at the beginning of the interval
is uniformly distributed between the two states
(i.e. :math:`\frac{1}{2}` prior probability for each state).

First, we want to compute the 

.. literalinclude:: in00.json
   :language: json
   :linenos:

The output reports the logarithm of the probability of jointly observing
state 0 at the beginning of the interval and
state 1 at the end of the interval.

.. literalinclude:: out00.json
   :language: json
   :linenos:

