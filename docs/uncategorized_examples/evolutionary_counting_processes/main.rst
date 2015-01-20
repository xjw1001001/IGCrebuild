table 1 by Minin and Suchard
============================

This example reproduces Table 1 of
"Counting labeled transitions in continuous-time
Markov models of evolution"
by Minin and Suchard, up to sampling error.

In the paper they sampled nucleotide sequences of length 1000
and computed labeled transition count expectations conditional
on those sequences, but here I've computed the expectations exactly
as if an infinite number of sequences had been sampled.
The discrepancies in the reported conditional expectations
seem to be within the sampling error.

.. literalinclude:: main.py
   :language: python
   :linenos:

.. literalinclude:: out.txt
   :linenos:
