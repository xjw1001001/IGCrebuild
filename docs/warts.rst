warts
=====

Despite its clear awesomeness, 
the jsonctmctree interface and its
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
as a todo-list, a bug list, or a list of 'issues.'

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
