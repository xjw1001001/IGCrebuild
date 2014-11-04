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

    states = j_in.get('root_posterior_states', None)
    expect = j_in.get('root_posterior_expect', None)
    none_count = sum(1 for x in (states, expect) if x is None)
    if none_count not in (0, 2):
        raise SimpleError('expected neither or both of '
                'root_posterior_states and '
                'root_posterior_expect to be provided')
    if not none_count:
        return np.array(states), np.array(expect, dtype=float)
    else:
        return None, None


def get_dwell_info(j_in):
    """These inputs are optional.
    """
    states = j_in.get('dwell_states', None)
    expect = j_in.get('dwell_expect', None)

The top-level input json object has the following members:
 * `node_count` : The number of nodes in the branching timeline.
   These include the root node, branching points, terminal points,
   and any additional points where observations are available along an edge.
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
 * `prior_feasible_states` : An array of prior feasible multivariate states
   at the root of the branching timeline.
   In our case, the root corresponds to a duplication event,
   so at the time of duplication the two variables in the bivariate
   model must share the same state.
   Therefore only four of the sixteen possible joint states are feasible.
 * `prior_distribution` : An array giving a prior probability for each
   of the listed feasible states.
 * `root_posterior_states` : Optional array requesting posterior
   linear combinations of probabilities of these multivariate states
   at the root node.
   The same linear combination will be used across all sites.
 * `root_posterior_expect` : Optional array of coefficients
   defining the linear combination of the multivariate states
   listed in `root_posterior_states`.
 * `dwell_states` : Optional array requesting posterior
   linear combinations of proportions of time spent in each
   of these multivariate states.
   The posterior linear combination will be reported
   for each site and each edge.
 * `dwell_expect` : Optional array of coefficients
   defining the linear combination of proportions of time
   spent in the the multivariate states listed in `root_posterior_states`.
 * `tree` : An object defining the shape of the branching timeline,
   the index of the process acting along each edge,
   and the edge-specific scaling factors of the process rates.
   The names `row` and `col` are in analogy to the representation
   of the branching timeline as a sparse matrix.
    * `row` : An array specifying the endpoint node closer to the root
      for each directed edge of the branching timeline.
    * `col` : An array specifying the endpoint node farther from the root
      for each directed edge of the branching timeline.
    * `rate` : An array specifying the rate scaling factor
      for each directed edge of the branching timeline.
    * `process` : An array specifying the index of the process
      acting on each directed edge of the branching timeline.
 * `processes` : An array of objects reprsenting processes
   that may act along one or more pre-specified edges of the branching timeline.
   Each object in the array is defined as follows:
    * `row` : An array of initial multivariate states for each
      feasible multivariate state transition.
      Because each multivariate state is defined by an array,
      this is an array of arrays.
    * `col` : An array of final multivariate states for each
      feasible multivariate state transition.
      Because each multivariate state is defined by an array,
      this is an array of arrays.
    * `rate` : An array of rates associated with each feasible
      multivariate state transition.
    * `expect` : An array of coefficients defining a linear combination
      of expected transition counts weighted according to the transition type.
      The output will include a value of the weighted sum
      for each edge at each site.
 * `observable_nodes` : An array of observable nodes.
   If more than one variable of the multivariate state is observable
   at the node, then the node index may be repeated in this list.
 * `observable_axes` : An array of observable axes.
   Together with the observable nodes array,
   these two arrays define which variables of the multivariate
   state are observable at which nodes.
   In our example, neither of the two variables of the bivariate process
   are observable at nodes that do not correspond to terminal nodes
   branching timeline, and at the terminal node that corresponds
   to the out-group, only one of the variables is observable.
 * `iid_observations` : An array of joint observations
   of variables at nodes on the branching timeline.
   Each observation is an array defining the variable state
   for each observable axis on each observable node.
   In our example, the number of observations corresponds
   to the length of the sequence alignment,
   and each observation has length 9 because we have
   four terminal nodes for which both variables are observable
   and one terminal node for which a single variable is observable.
