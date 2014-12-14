.. _example_9:

example 9
=========

An example of a molecular model of evolution of gene duplicates
using the following molecular data from primates.

.. literalinclude:: paralogs.fasta
   :linenos:

The maximum likelihood estimation uses
`L-BFGS-B
<http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>`_
with finite-differences approximations of the gradient of the log likelihood.

.. literalinclude:: mle.py
   :language: python
   :linenos:

.. literalinclude:: out00.txt
   :linenos:
