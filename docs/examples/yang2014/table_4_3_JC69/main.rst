.. _table_4_3_JC69:

table 4.3 JC69
==============

Although incomplete ambiguity of observations is not supported,
it is possible to emulate these types of observations using
an expanded multivariate state space.

In this example, we expand a univariate 4-state nucleotide process to a
:math:`2 \times 2 \times 2 \times 2` multivariate process.
The edge rate scaling factors are estimated using expectation maximization
to find the log likelihood of -6,262.01 reported in Ziheng Yang's 2014 textbook.

The output includes the results of 8 iterations of EM.


.. literalinclude:: em.py
   :language: python
   :linenos:

.. literalinclude:: em-out.txt
   :linenos:
