#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YDR438W --paralog2 YML018C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR438W --paralog2 YML018C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR438W --paralog2 YML018C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YDR438W --paralog2 YML018C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
