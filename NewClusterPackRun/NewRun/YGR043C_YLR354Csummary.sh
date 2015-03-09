#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YGR043C --paralog2 YLR354C --force False --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YGR043C --paralog2 YLR354C --force False --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YGR043C --paralog2 YLR354C --force True --clock True --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YGR043C --paralog2 YLR354C --force True --clock False --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
