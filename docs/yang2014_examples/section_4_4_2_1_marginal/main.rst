.. _section_4_4_2_1_marginal:

section 4.4.2.1 marginal reconstruction
=======================================

This is an implementation of the marginal reconstruction in section 4.4.2.1 of
`Molecular Evolution: A Statistical Approach`__, pages 127--128.

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

The output is probably transposed vs. the order you expect;
it gives an array of four arrays, each array representing a state,
whereas you might expect it to give nine arrays (one for each node)
with each of these nine arrays having length four.
The reason it is transposed this way is for indexing consistency,
which may have been a bad design choice.

In any case, you can see that the ancestral marginal reconstructions
in the json output agree with those given in the text::

    Node 0 (root):
        T : 0.055
        C : 0.901
        A : 0.037
        G : 0.007

    Node 6:
        T : 0.093
        C : 0.829
        A : 0.070
        G : 0.007

    Node 7:
        T : 0.153
        C : 0.817
        A : 0.026
        G : 0.004

    Node 8:
        T : 0.010
        C : 0.985
        A : 0.004
        G : 0.001
