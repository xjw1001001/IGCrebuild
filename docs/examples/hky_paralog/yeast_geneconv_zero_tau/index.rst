yeast HKY+geneconv with tau=0
=============================

In this example we estimate parameters of a molecular model
of the evolution of the DNA sequences of some paralogous genes
using the following molecular data from yeast.

.. literalinclude:: newick.tree
   :linenos:

.. literalinclude:: YDR502C_YLR180W.dat
   :linenos:

The molecular model is an HKY85 nucleotide substitution process
with the additional constraint that paralogous branches of the gene tree
have identical lengths.
There is no molecular clock constraint.


EM-like initialization plus quasi-Newton
----------------------------------------

The maximum likelihood search has two phases.
The first phase is an ad-hoc iterative procedure that tries to find
reasonable guesses of the parameter values,
and the second phase is a quasi-Newton search that uses
`L-BFGS-B
<http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>`_.
In this second phase, the derivatives of the log likelihood with respect
to edge-specific rate scaling factors are supplied using a closed-form
calculation, while the derivatives with respect to the mutational nucleotide
frequency parameters and the kappa parameter are supplied using
finite-difference approximations.

.. literalinclude:: main.py
   :language: python
   :linenos:

.. literalinclude:: main.sh
   :language: bash
   :linenos:

.. literalinclude:: out.txt
   :linenos:


EM for edge lengths only
------------------------

.. literalinclude:: inference02.py
   :language: python
   :linenos:

.. literalinclude:: out02.txt
   :linenos:
