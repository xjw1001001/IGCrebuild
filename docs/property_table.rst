property table
==============

This table indicates which reductions are available for which properties.

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


The core properties have the following meanings:

    LOGL
        log likelihood

    DERI
        derivatives of edge rate log scaling factors

    DWEL
        state occupancy distributions along edges

    TRAN
        transition count expectations

    ROOT
        root state distribution

    NODE
        state distributions at nodes


Reductions may be available for observations, edges, and states.
The reduction codes have the following meanings:

    D
        (D)istinct calculations -- no reduction

    S
        unweighted (S)ummation

    W
        (W)eighted summation

    N
        (N)ot applicable
