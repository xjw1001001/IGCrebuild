#!/bin/bash
python GenerateIndividualgBGCSummary.py --paralog1 YDR418W --paralog2 YEL054C --dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YDR418W --paralog2 YEL054C --dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YDR418W --paralog2 YEL054C --no-dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YDR418W --paralog2 YEL054C --no-dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
