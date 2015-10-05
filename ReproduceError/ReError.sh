#!/bin/bash
sbatch -o ReError-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/YLR406C_YDL075W_ReError.sh