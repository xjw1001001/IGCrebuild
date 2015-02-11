.. _example_8:

derivatives
===========

This example demonstrates the log likelihood derivatives request.
In particular, the request is for the derivatives of the
*log* likelihood with respect to the
*log* of the edge specific rate scaling factors.


2-state pure-death model with rate 5
------------------------------------

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


pure death with rate 5 + 1e-8
-----------------------------

Here we can see the same equality if we add a small 
delta to the rate scaling factor.

.. literalinclude:: in01.json
   :language: json
   :linenos:

.. literalinclude:: out01.json
   :language: json
   :linenos:


poisson with rate 3
-------------------

.. math::

    p \left( t \right)
    &= \frac{1}{2} \left( 1 + e^{-2 t} \right) \\
    \text{log} \left( p \left( 3 \right) \right)
    &\approx -0.690671495 \\
    \dfrac{
            \text{d}
        }
        {
            \text{d} \text{log}\left( t \right)
        }
    p \left( t \right)
    &= -\dfrac{2 t}{1 + e^{2 t}} \\
    \left.
    \dfrac{
            \text{d}
        }
        {
            \text{d} \text{log}\left( t \right)
        }
    p \left( t \right)
    \right|_{t=3}
    &= -\dfrac{6}{1 + e^{6}} \\
    &\approx -0.0148357


.. literalinclude:: in02.json
   :language: json
   :linenos:

.. literalinclude:: out02.json
   :language: json
   :linenos:


poisson with rate 3 + 1e-8
--------------------------

.. literalinclude:: in03.json
   :language: json
   :linenos:

.. literalinclude:: out03.json
   :language: json
   :linenos:

Using finite differences you can compute the derivative
of the log liklihood with respect to the scaling factor.
Three times this estimated derivative is approximately equal
to the computed derivative with respect to the log of the scaling factor.

.. math::

    3 * \dfrac{(-0.6906714954716674) - (-0.690671495422215)}{10^{-8}}
    &\approx 0.01483572

This follows from the chain rule of derivatives.
So if you want to use the interface to compute derivatives
of log likelihood with respect to edge rate scaling factors
instead of with respect to logs of edge rate scaling factors,
just divide the returned derivatives by the edge rate scaling factors.
