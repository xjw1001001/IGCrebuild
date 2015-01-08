.. _table_4_3_K80:

table 4.3 K80
=============

For a K80 model of molecular evolution along a fixed evolutionary tree shape,
use EM for maximum likelihood estimation of branch lengths and of the K80
kappa parameter given an observed alignment with IUPAC nucleotide ambiguity.

The output includes the results of 8 iterations of EM,
enough to converge to the maximum log likelihood -6,113.86 and the
maximum likelihood kappa parameter estiamate 3.561 given in
Table 4.3 of Ziheng Yang's 2014 textbook.


.. literalinclude:: main.py
   :language: python
   :linenos:

.. literalinclude:: out.txt
   :linenos:
