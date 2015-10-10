#!/bin/bash
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model MG94 --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model MG94 --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model MG94 --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model MG94 --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YJR048W --paralog2 YEL039C --model MG94 --no-clock --no-force --dir --gBGC --switch
