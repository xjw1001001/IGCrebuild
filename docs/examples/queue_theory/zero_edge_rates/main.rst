vanishing edge rates
====================

Edge rate scaling factors that are zero should not cause problems in general,
even when derivatives are computed.
In this example one of the edge rates is zero at an internal branch,
and two of the calculated derivatives are zero.


input
-----

.. literalinclude:: in00.json
   :language: json
   :linenos:

output
------

.. literalinclude:: out00.json
   :language: json
   :linenos:
