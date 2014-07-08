#!/bin/bash
sbatch -o CodonGeneconv-%j.out -p bigmem -w node93 test_run_script.sh