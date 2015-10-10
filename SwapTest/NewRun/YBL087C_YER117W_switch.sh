#!/bin/bash
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model HKY --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model MG94 --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model HKY --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model MG94 --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model HKY --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model MG94 --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model HKY --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model MG94 --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model HKY --no-clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YBL087C --paralog2 YER117W --model MG94 --no-clock --no-force --dir --gBGC --switch
