.. _example_6:

example 6
=========

This is from section 4.2.4 of
`Molecular Evolution: A Statistical Approach`__
by Ziheng Yang, 2014, Oxford University Press.

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


