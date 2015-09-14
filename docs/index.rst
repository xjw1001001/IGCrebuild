jsonctmctree overview
=====================

The purpose of this project is to provide a relatively easy way to query
properties of a continuous-time multivariate finite-state
Markov chain on a branching timeline with incomplete observations.

Some of these properties could include
the log likelihood,
its derivative with respect to logs of edge-specific rate scaling factors,
the expected proportions of time spent in subsets of states
along certain edges,
linear combinations of transition count expectations,
and posterior state distributions at nodes.

The properties themselves may be of interest,
or they may be only ingredients in a larger algorithm
that is out of scope for the jsonctmctree interface,
for example point estimation of parameter values.

Until a better user guide is written,
the best way to get started is by looking at one of the few examples below.


.. _reference:

reference
---------

.. toctree::
   :maxdepth: 2

   scene.rst
   reductions.rst
   property_table.rst
   examples/index.rst
   warts.rst
   about_the_docs.rst
   properties.rst
