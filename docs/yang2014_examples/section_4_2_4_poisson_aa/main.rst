.. _section_4_2_4_poisson_aa:

section 4.2.4 poisson amino acid model
======================================

This is from section 4.2.4 of
`Molecular Evolution: A Statistical Approach`__,
at the top of page 110.

__ http://abacus.gene.ucl.ac.uk/MESA/

It computes the log likelihood for a poisson model of molecular evolution
for an evolutionary tree of primates.

First, here is the json input that provides
a separate observation for each column in the alignment
and computes the log likelihood.

.. literalinclude:: in00.json
   :language: json
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :linenos:


Here is the json input that
reduces the input by considering only multiplicities of distinct patterns
of equality or inequality among amino acids at alignment columns.

.. literalinclude:: in01.json
   :language: json
   :linenos:

.. literalinclude:: out01.json
   :language: json
   :linenos:

The conditional expected number of transitions per edge 
can be used to iteratively update the parameters corresponding to
edge-specific rate scaling factors,
reaching the -16,566.60 log likelihood after four iterations.
This is an expectation maximization.

.. literalinclude:: em.py
   :language: python
   :linenos:

.. literalinclude:: em_out.json
   :language: json
   :linenos:

And using the 'extras' module:

.. literalinclude:: extras_em.py
   :language: python
   :linenos:

.. literalinclude:: extras_em.out


