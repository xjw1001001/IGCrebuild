.. _section_4_2_4_mtmam:

section 4.2.4 MTMAM
===================

This is another example from section 4.2.4 of
`Molecular Evolution: A Statistical Approach`__, page 109.

__ http://abacus.gene.ucl.ac.uk/MESA/

It computes the log likelihood for the empirical MTMAM
model of molecular evolution for an evolutionary tree of primates.

The conditional expected number of transitions per edge 
divided by the conditional expected exit rate from each state
can be used to iteratively update the parameters corresponding to
edge-specific rate scaling factors,
reaching the -14,558.59 log likelihood after six iterations.
This is an expectation maximization.

.. literalinclude:: em.py
   :language: python
   :linenos:

.. literalinclude:: em_out.json
   :language: json
   :linenos:

And using the generic edge rate EM in the 'extras' module:

.. literalinclude:: em_extras.py
   :language: python
   :linenos:

.. literalinclude:: em_extras.out
