example 3
=========

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

input
^^^^^

.. literalinclude:: in00.json
   :language: json
   :linenos:


output
^^^^^^

.. literalinclude:: out00.json
   :language: json
   :linenos:


estimating parameters
---------------------

EM iteration 1 input
^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in01.json
   :language: json
   :linenos:


EM iteration 1 output
^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: out01.json
   :language: json
   :linenos:


EM iteration 2 input
^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: in02.json
   :language: json
   :linenos:


EM iteration 2 output
^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: out02.json
   :language: json
   :linenos:
