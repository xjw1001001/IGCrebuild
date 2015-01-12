.. _example_4_1:

example 4.1
===========

This is example 4.1 in
`Molecular Evolution: A Statistical Approach`__, section 4.2.2.1, page 105.

__ http://abacus.gene.ucl.ac.uk/MESA/

The stochastic process is called K80_
whose state space is {T, C, A, G}
which we will encode as {0, 1, 2, 3}.
The example uses :math:`\kappa = 2`,
but with a slightly different parameterization of the rate matrix,
so all rates in this example are scaled by 0.25
relative to the parameterization on Wikipedia.

.. _K80: http://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_.28Kimura.2C_1980.29.5B2.5D

.. literalinclude:: in00.json
   :language: json
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :linenos:
