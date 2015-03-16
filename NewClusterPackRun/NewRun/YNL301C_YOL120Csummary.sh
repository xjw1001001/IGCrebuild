#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YNL301C --paralog2 YOL120C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL301C --paralog2 YOL120C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL301C --paralog2 YOL120C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YNL301C --paralog2 YOL120C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
