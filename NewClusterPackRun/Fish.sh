#!/bin/bash
sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/TLR5a_TLR5b_switch.sh