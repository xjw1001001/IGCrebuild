#!/bin/bash
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model HKY --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model HKY --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model HKY --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model HKY --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model HKY --no-clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --dir --gBGC --switch
