#!/bin/bash
python GenerateIndividualDirSummary.py --paralog1 YLR333C --paralog2 YGR027C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR333C --paralog2 YGR027C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR333C --paralog2 YGR027C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR333C --paralog2 YGR027C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
