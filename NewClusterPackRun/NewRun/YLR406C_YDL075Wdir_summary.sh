#!/bin/bash
python GenerateIndividualDirSummary.py --paralog1 YLR406C --paralog2 YDL075W --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR406C --paralog2 YDL075W --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR406C --paralog2 YDL075W --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YLR406C --paralog2 YDL075W --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
