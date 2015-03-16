#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YBL087C --paralog2 YER117W --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
