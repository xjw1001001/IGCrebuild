#!/bin/bash
python GenerateIndividualgBGCSummary.py --paralog1 YNL049C --paralog2 YIL109C --dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL049C --paralog2 YIL109C --dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL049C --paralog2 YIL109C --no-dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL049C --paralog2 YIL109C --no-dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
