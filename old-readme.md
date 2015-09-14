Introduction
============

The technical motivation of this project was to work with likelihoods
for parametric models of multivariate finite-state continuous-time evolution
on a fixed branching history with branch-specific rate scaling factors
that are also parameters in the model.
The evolution along different branches may be governed by different
predetermined processes, as long as all such processes
share the same state space.
At each node of the branching process,
each axis the multivariate state space may be either completely
observed or completely unobserved.

The rate matrices describing the continuous-time evolution
of the multivariate state are assumed to be medium-sized and sparse,
with a state space defined by the cartesian product of a small number of
small finite discrete subspaces.

For this project, the size of the joint state space should be
somewhere between hundreds of states and tens of thousands of states.
For smaller state spaces, existing software packages would be more appropriate.
This project is also not applicable to models with state spaces
that are cartesian products of a large number of large subspaces
resulting in a combinatorial explosion of the state space size.


Examples from molecular evolution, with categorical variables and branching
===========================================================================


Muse-Gaut 1994 (MG94) codon model of molecular sequence evolution
----------------------------------------------------------

This model considers the joint stochastic evolution of the three
nucleotide positions in a codon.
In this model each of the three nucleotide positions evolves according to
the Hasegawa-Kishino-Yano 1985 (HKY85) model,
but with a couple of dependence mechanisms.
The first factor introducing dependence is to disallow
the three variables from having states that jointly code for a stop codon.
The second factor is that rates of nucleotide change that affect the amino
acid encoded by joint states are scaled differently
from nucleotide changes that do not change the encoded amino acid.

The multivariate sample space is {A, C, G, T}^3 although multivariate
states corresponding to stop codons are effectively excluded.


Muse 1995 di-nucleotide model of RNA pairing constrained DNA evolution
----------------------------------------------------------------------

This bivariate generalization HKY85 models the co-evolution
of nucleotides at positions that are nominally paired according to RNA
secondary structure.
A parameter controlling the tendency to favor paired states introduces
the dependence between variables.

The multivariate sample space is {A, C, G, T}^2.


Covarion-like or Tuffley-Steel models of molecular evolution
------------------------------------------------------------

These kinds of models track an extra latent binary variable
together with the primary observable variable.
One example is an evolving latent binary variable that controls
the overall rate of evolution at a nucleotide position
that evolves according to HKY85.

The multivariate sample space is {0, 1} x {A, C, G, T}.


Examples from other fields, with bounded integer variables and no branching
===========================================================================

The following examples are bivariate stochastic models
on continuous time intervals.
The sample space of each variable is a finite range of integers
which is interpreted as a range of counts rather than as a set of categories,
although this distinction is not used in the analysis for this project.
Other approaches that have the goal of answering similar
questions may care about this distinction,
for example approaches that use mean field theory or diffusion approximations.

Unlike the molecular evolution models that assume a branching timeline
with observations at the tips of the tree,
these models assume a linear timeline with observations
at arbitrary points along the interval.
As in the models of molecular evolution,
it is possible that only a subset of the variables
in the multivariate process are observable.


That economics model with finite discrete fixed size inbox and outbox
---------------------------------------------------------------------

You have an in-box and an out-box each with a capacity to hold
a finite number of parts and widgets respectively.
The in-box receives parts stochastically at some expected rate,
and the out-box delivers widgets stochastically at some other expected rate.
A worker creates widgets from parts stochastically at some rate
if parts are available in the in-box
and if the out-box has room for the widget.
The rate parameters are assumed to be constant over time.

The multivariate sample space is {0, 1, ..., n} x {0, 1, ..., k}
where n and k are the inbox and outbox sizes respectively.


Finite discrete compartmental models in epidemiology
----------------------------------------------------

Say that some number of people are infected with some disease of interest
within a total population of fixed size, and that some other number
of people within that population are recovered from the disease.

We assume that these numbers change stochastically in discrete steps
over continuous time, according to dynamics defined by rate parameters
that are constant over time.

The multivariate sample space is {0, 1, ..., n} x {0, 1, ..., n}
where n is the total population, although the joint states
that corresponds to a sum of infected and recovered that is larger
that the total population are effectively excluded.


Stochastic treatment of finite discrete Predator-Prey models in ecology
-----------------------------------------------------------------------

We assume that the stochastic increase and decrease of the predator and prey
populations occurs in discrete steps over continuous time
according to expected rates that are time-homogeneous functions
of the two finite population sizes.

The multivariate sample space is {0, 1, ..., n} x {0, 1, ..., k}
where n and k are the two maximum population sizes.
