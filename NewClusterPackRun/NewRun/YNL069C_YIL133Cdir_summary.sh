#!/bin/bash
python GenerateIndividualDirSummary.py --paralog1 YNL069C --paralog2 YIL133C --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL069C --paralog2 YIL133C --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL069C --paralog2 YIL133C --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualDirSummary.py --paralog1 YNL069C --paralog2 YIL133C --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
