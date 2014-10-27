In this documentation file we explain an example input file.

The input file is in the [JSON](http://json.org/) meta-format.
It can be read and written by hand and by standard library functions
available for most programming languages.
As a meta-format, JSON has a similar but not identical niche as
[XML](http://en.wikipedia.org/wiki/XML) and
[RDF](http://en.wikipedia.org/wiki/Resource_Description_Framework).

JSON consists of recursively nested objects and arrays,
with terminal values consisting of strings, numbers, true, false, and nil.

The input consists primarily of information about the structure
of the model, specific rates associated with state transitions
and with edges of the branching timeline, and observations
at points on the branching timeline.

The top-level input json object has the following members:
 * `node_count` : The number of nodes in the branching timeline.
   These can be branching points, or terminal points,
   or points where observations are available along an edge.
   In our example, there are nine nodes.
 * `process_count` : The number of different processes.
   Each edge between nodes on the branching timeline is assumed
   to evolve according to exactly one of these processes.
   In our example, there are two distinct processes, corresponding to
   presence or absence of a gene duplicate.
 * `state_space_shape` : This is an array of integers defining
   the shape of the multivariate state space.  Our example model
   is bivariate with each variable having four possible states,
   so the multivariate state space shape is `[4, 4]`.
 * TODO: document more inputs...
