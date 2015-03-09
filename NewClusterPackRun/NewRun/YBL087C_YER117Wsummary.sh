#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force False --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force False --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force True --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force True --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
