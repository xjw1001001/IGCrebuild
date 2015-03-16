#!/bin/bash
python GenerateIndividualDirSummary.py --paralog1 YNL049C --paralog2 YIL109C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL049C --paralog2 YIL109C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL049C --paralog2 YIL109C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL049C --paralog2 YIL109C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
