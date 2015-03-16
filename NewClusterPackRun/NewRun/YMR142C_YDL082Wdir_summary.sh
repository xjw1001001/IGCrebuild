#!/bin/bash
python GenerateIndividualDirSummary.py --paralog1 YMR142C --paralog2 YDL082W --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YMR142C --paralog2 YDL082W --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YMR142C --paralog2 YDL082W --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YMR142C --paralog2 YDL082W --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
