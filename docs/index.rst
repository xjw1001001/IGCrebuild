jsonctmctree overview
=====================

The purpose of this project is to provide a relatively easy way to query
properties of a continuous-time multivariate finite-state
Markov chain on a branching timeline with incomplete observations.

The basic idea is to pass a
:ref:`scene`
and some
:ref:`requests`
like

.. code-block:: json

    {
        "scene" : {},
        "requests" : []
    }

to a library function, command-line program, or website
that implements the jsonctmctree interface,
and it will respond with something like

.. code-block:: json

    {
        "status" : "",
        "reponses" : []
    }



.. _examples:

examples
--------

.. toctree::
   :maxdepth: 2

   tut00/main.rst
   tut01/main.rst
   tut02/main.rst
   tut03/main.rst


.. _reference:

reference
---------

.. toctree::
   :maxdepth: 2

   scene.rst
   reductions.rst
   property_table.rst
   about_the_docs.rst


.. _extended_properties:

extended properties
-------------------

.. toctree::
   :maxdepth: 2

   properties.rst
