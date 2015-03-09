#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YPR159W --paralog2 YGR143W --force False --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YPR159W --paralog2 YGR143W --force False --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YPR159W --paralog2 YGR143W --force True --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YPR159W --paralog2 YGR143W --force True --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
