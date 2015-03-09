#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YNL069C --paralog2 YIL133C --force False --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL069C --paralog2 YIL133C --force False --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL069C --paralog2 YIL133C --force True --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL069C --paralog2 YIL133C --force True --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
