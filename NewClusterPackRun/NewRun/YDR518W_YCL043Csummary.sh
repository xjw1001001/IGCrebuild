#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YDR518W --paralog2 YCL043C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR518W --paralog2 YCL043C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR518W --paralog2 YCL043C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR518W --paralog2 YCL043C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
