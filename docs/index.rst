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

Here's a cryptic table for you to puzzle over.

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


user documentation
------------------

.. toctree::
   :maxdepth: 2

   tut00/main.rst
   tutorial_1.rst
   about_the_docs.rst
