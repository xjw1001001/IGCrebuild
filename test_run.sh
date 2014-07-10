#!/bin/bash
sbatch -o CodonGeneconv-%j.out -p bigmem -w node93 test_run_script.sh
sbatch --mail-type=END --mail-type=FAIL --mail-user=xji3@ncsu.edu wait.bash