LOGL
====

Log likelihoods.


DNNLOGL
-------

    request
        This property request applies no reduction.

    response
        The response is a 1d array
        of observation-specific log likelihoods.
        The length of the array is equal to the
        number of iid observations provided in the scene.


SNNLOGL
-------

    request
        This property request applies a pre-defined reduction,
        summing the log likelihoods over all iid observations.

    response
        The response is a single number,
        representing the sum of log likelihoods over iid observations.


WNNLOGL
-------

    request
        This property request applies a user-defined observation reduction,
        computing a weighted sum of log likelihoods over iid observations.

    response
        The response is a single number,
        representing a weighted sum of log likelihoods over iid observations.
