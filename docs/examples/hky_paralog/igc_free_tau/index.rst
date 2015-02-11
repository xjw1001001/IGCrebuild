HKY+IGC free tau
================

An example of a molecular model of evolution of gene duplicates
using the following molecular data from primates.

.. literalinclude:: paralogs.fasta
   :linenos:


finite differences gradient
---------------------------

The following maximum likelihood estimation uses
`L-BFGS-B
<http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>`_
with finite-differences approximations of the gradient of the log likelihood.

.. literalinclude:: mle.py
   :language: python
   :linenos:

.. literalinclude:: out00.txt
   :linenos:


explicit edge rate derivatives
------------------------------

This more complicated implementation still uses
`L-BFGS-B
<http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>`_
but with the edge-rate component of the gradient computed
using calculus instead of finite differences.
The entries of the gradient corresponding to the other parameters
of the model (mutational nucleotide distribution, transition/transversion
rate ratio, and gene conversion rate) are still estimated
using a finite differences approximation.

.. literalinclude:: mle-gradient.py
   :language: python
   :linenos:

.. literalinclude:: out01.txt
   :linenos:
