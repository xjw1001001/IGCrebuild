"""
This is a helper script.

The idea is to use codeml to find maximum likelihood estimates
for a Muse-Gaut 1994 codon model with a tree of paralogs.
We cannot constrain codeml to use the same branch lengths
for the corresponding branches of the two paralog trees,
so instead we will just average the maximum likelihood lengths
of the corresponding branches.

The codeml analysis itself takes place within a temporary directory.
Files are copied to that directory and read from that directory.

"""
from __future__ import print_function, division

import subprocess
import contextlib
import argparse
import tempfile
import os
import sys
import shutil

import numpy as np
from numpy.testing import assert_equal

import dendropy

_codeml_ctl_template = """
      seqfile = {alignment}
     treefile = {tree}
      outfile = mlc

        noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1
    * CodonFreq = 0
      * estFreq = 1

    CodonFreq = 4
      estFreq = 1

        ndata = 1
        clock = 0
       aaDist = 0

        model = 0

      NSsites = 0 

        icode = 0
        Mgene = 0

    fix_kappa = 0
        kappa = 5.86954
    fix_omega = 0
        omega = 0.087136

    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 5

        getSE = 1
 RateAncestor = 0

   Small_Diff = 5e-7
    cleandata = 0
  fix_blength = 0
       method = 0

"""


# This is an example of the codeml output
# which by default goes into a file named mlc.
"""

...

tree length =   1.65304

(((((((6: 0.034271, 8: 0.032884): 0.007791, 7: 0.168901): 0.000004, 4: 0.092079): 0.012583, 10: 0.071975): 0.100903, 3: 0.218238): 0.000004, (((((5: 0.081158, 12: 0.025724): 0.015747, 11: 0.072845): 0.000004, 9: 0.123571): 0.032967, 13: 0.065981): 0.096087, 2: 0.157856): 0.000004): 0.119675, 1: 0.121793);

(((((((cerevisiaeYDR450W: 0.034271, paradoxusYDR450W: 0.032884): 0.007791, mikataeYDR450W: 0.168901): 0.000004, kudriavzeviiYDR450W: 0.092079): 0.012583, bayanusYDR450W: 0.071975): 0.100903, castelliiYDR450W: 0.218238): 0.000004, (((((cerevisiaeYML026C: 0.081158, paradoxusYML026C: 0.025724): 0.015747, mikataeYML026C: 0.072845): 0.000004, kudriavzeviiYML026C: 0.123571): 0.032967, bayanusYML026C: 0.065981): 0.096087, castelliiYML026C: 0.157856): 0.000004): 0.119675, kluyveriYML026C: 0.121793);

Detailed output identifying parameters

kappa (ts/tv) =  4.87425

Frequency parameters:
       0.30051 (T)   0.20980 (C)   0.29395 (A)   0.19574 (G)

...

omega (dN/dS) =  0.05900

...

"""


def get_prefixed_line_index(lines, prefix):
    indices = []
    for i, line in enumerate(lines):
        if line.startswith(prefix):
            indices.append(i)
    if len(indices) == 1:
        return indices[0]
    else:
        raise Exception('found %d lines beginning with the prefix "%s"' % (
            (len(indices), prefix)))


def scrape_codeml_output(lines):

    # Look for indices of lines that begin with certain prefixes.
    tree_length_idx = get_prefixed_line_index(lines, 'tree length =')
    omega_idx = get_prefixed_line_index(lines, 'omega (dN/dS) =')
    kappa_idx = get_prefixed_line_index(lines, 'kappa (ts/tv) =')
    nt_freq_idx = get_prefixed_line_index(lines, 'Frequency parameters:')

    # Extract the simplest parameter values.
    # The tree length is not really a parameter but we can get it anyway.
    kappa = float(lines[kappa_idx].split()[-1])
    omega = float(lines[omega_idx].split()[-1])
    tree_length = float(lines[tree_length_idx].split()[-1])

    # Extract the mutational nucleotide distribution parameters.
    # 0.30051 (T)   0.20980 (C)   0.29395 (A)   0.19574 (G)
    t, tx, c, cx, a, ax, g, gx = lines[nt_freq_idx + 1].split()
    assert_equal((tx, cx, ax, gx), ('(T)', '(C)', '(A)', '(G)'))
    pi = np.array([float(x) for x in (a, c, g, t)])

    # Extract the newick tree which also contains branch lengths.
    i = tree_length_idx
    integer_node_tree_idx = i + 2
    string_node_tree_idx = i + 4
    newick_tree = lines[string_node_tree_idx]

    # Return the values.
    return pi, kappa, omega, newick_tree


# Use this to temporarily change to a directory.
@contextlib.contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)


def main(args):

    # Define the contents of codeml.ctl.
    codeml_ctl_contents = _codeml_ctl_template.format(
            alignment=args.alignment,
            tree=args.tree)

    # Require that a file named codeml.ctl
    # does not exist in the current directory.
    if os.path.exists('codeml.ctl'):
        raise Exception('codeml.ctl already exists in the current directory')

    # Create a control file in the current directory.
    with open('codeml.ctl', 'w') as fout:
        fout.write(codeml_ctl_contents)

    # Report contents of the current directory.
    print('contents of the current directory:')
    print(os.listdir('.'))

    # Run codeml in the current directory.
    proc = subprocess.Popen(
            args=(args.codeml,),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    stdoutdata, stderrdata = proc.communicate()
    print('stdoutdata:')
    print(stdoutdata)
    print('stderrdata:')
    print(stderrdata)

    # Read the inferred parameter values from the output.
    with open('mlc') as fin:
        lines = fin.readlines()
    pi, kappa, omega, newick_tree = scrape_codeml_output(lines)

    # Print the inferred parameter value estimates.
    print('parameter estimates')
    print('-------------------')
    print('mutational nucleotide frequencies:')
    print(pi)
    print('kappa:', kappa)
    print('omega:', omega)
    print('newick tree with branch length estimates:')
    print(newick_tree)


if __name__ == '__main__':

    # Get the command-line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('--codeml', required=True,
            help='name of the codeml binary including the full path')
    parser.add_argument('--alignment', required=True,
            help=(
                'full path to a file named something like '
                'YLR284C_YOR180C_input.fasta'))
    parser.add_argument('--tree', required=True,
            help=(
                'full path to a newick file that has exactly one sequence '
                'of the alignment at each tip of the tree'))
    args = parser.parse_args()

    # Create a temporary directory,
    # and copy the alignment and tree into this directory.
    tmpdir = tempfile.mkdtemp()
    shutil.copy(args.alignment, tmpdir)
    shutil.copy(args.tree, tmpdir)

    # Run the script within a temporary directory.
    print('current directory before analysis:', file=sys.stderr)
    print(os.getcwd(), file=sys.stderr)
    with cd(tmpdir):
        print('current directory during analysis:', file=sys.stderr)
        print(os.getcwd(), file=sys.stderr)
        main(args)
    print('current directory after analysis:', file=sys.stderr)
    print(os.getcwd(), file=sys.stderr)
