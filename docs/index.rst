.. jsonctmctree documentation master file, created by
   sphinx-quickstart on Tue Dec  2 20:58:55 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

jsonctmctree overview
=====================

jsonctmctree will be among other things a specification
of a json interface for querying properties of a continuous-time
finite-multivariate-state Markov chain on a branching timeline
with incomplete observations.

+------+---------------+---------------+---------------+
|      |  observations |     edges     |    states     |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
|      | D | S | W | N | D | S | W | N | D | S | W | N | 
+======+===+===+===+===+===+===+===+===+===+===+===+===+
| LOGL | Y | Y | Y | . | . | . | . | Y | . | . | . | Y |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
| DERI | Y | Y | Y | . | Y | . | . | . | . | . | . | Y |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
| DWEL | Y | Y | Y | . | Y | . | Y | . | Y | . | Y | . |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
| TRAN | Y | Y | Y | . | Y | Y | . | . | . | . | . | Y |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
| ROOT | Y | Y | Y | . | . | . | . | Y | Y | . | Y | . |
+------+---+---+---+---+---+---+---+---+---+---+---+---+
| NODE | Y | Y | Y | . | . | . | . | Y | Y | . | Y | . |
+------+---+---+---+---+---+---+---+---+---+---+---+---+


Contents:

.. toctree::
   :maxdepth: 2

   tutorial_0.rst
   tutorial_1.rst
   about_the_docs.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

