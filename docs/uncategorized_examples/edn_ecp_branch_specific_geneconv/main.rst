edge-specific posterior estimates
=================================

Compute edge-specific posterior parameter estimates
using ratios of posterior expectations.

This uses the HKY85 model of evolution of gene duplicates
using the following molecular data from primates,
with the added twist that an effect that 'homogenizes'
paralogous sequences within a species is present.

.. literalinclude:: paralogs.fasta
   :linenos:

.. literalinclude:: main.py
   :language: python
   :linenos:

.. literalinclude:: out.txt
