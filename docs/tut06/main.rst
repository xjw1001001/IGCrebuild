.. _example_7:

example 7
=========

This is another example from section 4.2.4 of
`Molecular Evolution: A Statistical Approach`__
by Ziheng Yang, 2014, Oxford University Press.

__ http://abacus.gene.ucl.ac.uk/MESA/

It computes the log likelihood for a poisson model of molecular evolution
for an evolutionary tree of primates,
using the empirical MTMAM model.

The conditional expected number of transitions per edge 
divided by the conditional expected exit rate from each state
can be used to iteratively update the parameters corresponding to
edge-specific rate scaling factors,
reaching the -14,558.59 log likelihood after seven iterations.
This is an expectation maximization.

.. literalinclude:: em.py
   :language: python
   :linenos:

.. literalinclude:: em_out.json
   :language: json
   :linenos:


