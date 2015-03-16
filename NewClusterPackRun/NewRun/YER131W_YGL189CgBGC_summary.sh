#!/bin/bash
python GenerateIndividualgBGCSummary.py --paralog1 YER131W --paralog2 YGL189C --dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YER131W --paralog2 YGL189C --dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YER131W --paralog2 YGL189C --no-dir --no-force --clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
python GenerateIndividualgBGCSummary.py --paralog1 YER131W --paralog2 YGL189C --no-dir --no-force --no-clock --model MG94 --sump ./NewPackageNewRun/ --pairp ./NewPackageNewRun/
