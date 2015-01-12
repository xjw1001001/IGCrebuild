"""
Reproduce parameter estimates and log likelihoods
from an example in Ziheng Yang's 2014 textbook.

Section 4.7. Page 144.
plastid rbcL genes from 12 plant species.
1428 nucleotide sites.
A fixed tree shape is used.
The models are time-reversible so the root of the tree does not matter.

Table 4.3.
The number of degrees of freedom is indicated by p.
The log likelihood is indicated by l.
The model modifier +C indicates that the nucleotide columns
of the coding sequence alignment are partitioned
according to the offset {0, 1, 2} within the codon.

model: JC69
p: 21
l: -6262.01

model: K80
p: 22
l: -6113.86
k: 3.561

model: HKY85
p: 25
l: -6101.76
k: 3.620

model: JC69+C
p: 23
l: -5922.76
r1: 1
r2: 0.556
r3: 5.405

model: K80+C
p: 26
l: -5728.76
k1: 1.584
k2: 0.706
k3: 5.651
r1: 1
r2: 0.556
r3: 5.611

model: HKY85+C
p: 35
l: -5624.70
k1: 1.454
k2: 0.721
k3: 6.845
r1: 1
r2: 0.555
r3: 5.774

"""
