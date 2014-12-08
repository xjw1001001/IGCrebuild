weighted reductions
===================

To apply reductions that are more complicated than summation,
the coefficients (or 'weights') of the linear combination
need to be specified.


.. _observation_reduction:

observation_reduction
---------------------

Specifies a sparse weighted sum over observations.

    observation_indices : 1d array of integers
        Each observation index in the array
        points to a row of the iid_observations array in the scene object.

    weights : 1d array of numbers
        Coefficients of the custom linear combination.


.. _edge_reduction:

edge_reduction
--------------

Specifies a sparse weighted sum over edges.

    edges : 1d array of integers
        Each edge index in the array
        points to an edge of the branching timeline.

    weights : 1d array of numbers
        Coefficients of the custom linear combination.


.. _state_reduction:

state_reduction
---------------

Specifies a sparse weighted sum over multivariate states.

    states : 2d array of integers
        An array of multivariate states.

    weights : 1d array of numbers
        Coefficients of the custom linear combination.


.. _transition_reduction:

transition_reduction
--------------------

Specifies a sparse weighted sum over state transitions.

    row_states : 2d array of integers
        An array of multivariate states.

    column_states : 2d array of integers
        An array of multivariate states.

    weights : 1d array of numbers
        Coefficients of the custom linear combination.

