From the paml manual:

You can collect initial values into a file called [...] in.codeml if
you are running codeml.  This file should contain as many numbers,
separated by white spaces,
as the number of parameters that are being optimized by the program.
So if the program is estimating 56 parameters (say 51 branch lengths,
1 kappa, and 5 other parameters from the omega distribution),
you should put 56 numbers in the file.

The order is the same as in the main output (below the lnL line for each tree.

One way of generating the in.codeml file is to run a data set, and then copy
initial values from the rub file or from the main output file.
The rub file records the iteration process and has one line
for each round of iteration.
When you run the program, look at lnL0 printed out on the screen
and check that it is the same as recorded in rub.

If you have already obtained parameter estimates before and do not want the
program to re-estimate them and only want to do some analysis based
on those estiamtes such as reconstructing ancestral sequences,
insert -1 before the initial values.

In the file of initial values, if you use -1 at the start, the program
assumes the original variables, while if you don't, the program assumes
transformed variables.
