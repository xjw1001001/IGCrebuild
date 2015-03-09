#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YMR142C --paralog2 YDL082W --force False --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YMR142C --paralog2 YDL082W --force False --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YMR142C --paralog2 YDL082W --force True --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YMR142C --paralog2 YDL082W --force True --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
