#!/bin/bash
python GenerateIndividualgBGCSummary.py --paralog1 YNL069C --paralog2 YIL133C --dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL069C --paralog2 YIL133C --dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL069C --paralog2 YIL133C --no-dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YNL069C --paralog2 YIL133C --no-dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
