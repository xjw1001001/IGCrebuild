warts
=====

The jsonctmctree interface and its
reference implementation are not perfect.
Some of this awkwardness is caused by bad planning and bad design on one hand,
and on the other hand by the deliberate tradeoffs that seek to balance
the goals of maximizing expressive power, maximizing speed,
minimizing memory usage, minimizing the sizes of the inputs and outputs,
maximizing human readability of the inputs and outputs,
maximizing accessibility from multiple programming languages,
and generally minimizing complexity.

With that said, here's a compilation of 'warts.'
It's under construction and should not necessarily be understood
as a todo-list, a bug list, or a list of 'issues' or a FAQ.

* Error reporting is not so great, in general.
* Negative rates do not currently cause an informative error.
* Coefficients of linear combinations are named as 'weights'
  which implies non-negativity, but in principle there is no reason
  for such a non-negativity constraint.
* Zeros are not currently sparsified automatically.
* The shapes of the input arrays are unnecessarily tightly constrained.
  For example, the row_states and column_states in the definition
  of transition rate matrices are required to be 2d,
  but in principle a transition matrix devoid of rates (all rates zero)
  should be allowed, and its row_states and column_states should
  both be [] which is 1d with length 0.
* Observations are only allowed at points,
  not along intervals with positive length.
  An example of a currently unsupported but theoretically possible feature
  would be continuous observation of one variable on an interval,
  jointly evolving with an unobserved latent variable.
* Although aggregate observations are representable through the interface,
  they are currently not efficiently supported
  in the reference implementation.
  For example you might imagine simplifying your representation of a codon
  model by tacking on an extra co-evolving variable corresponding to the
  translated amino acid.
  This would theoretically explode the size of the state space to
  61*20 or 64*20 depending on the codon representation,
  but only 1/20 of this space would actually be used.
  Currently the reference implementation is not clever enough to
  avoid constructing arrays whose size is equal to that of the
  joint state space.
* Noisy observations are not supported by the interface.
  Flexible representation of noisy observations could enable many
  extra features, for example it could represent the sampling uncertainty
  given a hidden population state.
  Or it could represent an actual measurement uncertainty.
* Partitions of the observations are not explicitly supported
  by the interface, but they could be emulated in a few ways.
  One way could be through the observation reduction in each request.
  Another way would be to split the partitions across function calls.
* Second derivatives and cross derivatives
  of the log likelihood with respect to logs of edge rate scaling factors
  are not supported.
  This is theoretically possible to add and would be efficient to calculate,
  and would be more accurate than finite differences.
* The interface uses edge rate scaling factors,
  but derivatives of the log likelihood are taken with respect
  to the *logs* of the edge rate scaling factors.
  Presumably most users would not care about derivatives,
  and it would be more confusing to require that they provide
  logs of scaling factors instead of scaling factors.
  On the other hand, for complicated situations it is easier to work with logs
  of edge rate scaling factors, and if you are asking for derivatives
  then you are probably already dealing with a complicated situation,
  so the extra complication of this log transform would be the least
  of your worries.
* Uncertainty in the branching structure of the timeline
  is not supported.
* Random sampling is not supported.
* Models that include random effects are not supported.
* Inference and hypothesis testing are not supported,
  although they could easily be built on top of this interface.
  On the other hand this layering would result in algorithms
  that are slower than would be possible with more integrated solutions.
* Bayesian inference is not supported.
* Mixture models are not treated specially,
  but they could be emulated using block diagonal rate matrices.
* The extra structure of time-reversible rate matrices is not used.
