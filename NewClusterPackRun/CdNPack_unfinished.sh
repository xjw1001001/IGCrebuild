#!/bin/bash
sbatch -o summary-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/YMR143W_YDL083C_clock_gBGC_unfinished.sh
