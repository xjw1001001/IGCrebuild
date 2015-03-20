#!/bin/bash
sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/YIR033W_YKL020C_clock_dir_gBGC_unfinished.sh
sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/YIR033W_YKL020C_clock_gBGC_unfinished.sh
