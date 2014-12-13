.. _example_8:

example 8
=========

Examine the log likelihood derivatives request.
This may be confusing, because the derivative is of the 
*log* likelihood with respect to the
*log* of the edge specific rate scaling factor.


rate 5
------

In this example with a single branch,
the log likelihood is equal to
the derivative of the log likelihood with respect to
the log scaling factor.

.. math::

    p
    &= e^{-t} \\
    \text{log}\left( p \right)
    &= -t \\
    \dfrac{
            \text{d}
        }
        {
            \text{d} \text{log}\left( t \right)
        }
    \text{log} \left( p \right)
    &= -t


.. literalinclude:: in00.json
   :language: json
   :linenos:

.. literalinclude:: out00.json
   :language: json
   :linenos:


rate 5 + 1e-8
-------------

Here we can see the same equality if we add a small 
delta to the rate scaling factor.

.. literalinclude:: in01.json
   :language: json
   :linenos:

.. literalinclude:: out01.json
   :language: json
   :linenos:
