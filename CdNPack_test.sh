#!/bin/bash
sbatch -o cd-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./NewRun/YLR406C_YDL075Wcd.sh
