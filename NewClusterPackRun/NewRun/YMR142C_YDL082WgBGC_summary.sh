#!/bin/bash
python GenerateIndividualgBGCSummary.py --paralog1 YMR142C --paralog2 YDL082W --dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YMR142C --paralog2 YDL082W --dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YMR142C --paralog2 YDL082W --no-dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YMR142C --paralog2 YDL082W --no-dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
