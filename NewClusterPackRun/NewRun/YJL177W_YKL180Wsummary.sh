#!/bin/bash
python GenerateIndvidualSummary.py --paralog1 YJL177W --paralog2 YKL180W --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YJL177W --paralog2 YKL180W --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YJL177W --paralog2 YKL180W --force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndvidualSummary.py --paralog1 YJL177W --paralog2 YKL180W --force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
