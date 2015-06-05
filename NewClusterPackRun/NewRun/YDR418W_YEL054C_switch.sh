#!/bin/bash
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model MG94 --clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model MG94 --clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model MG94 --clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model MG94 --clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YDR418W --paralog2 YEL054C --model MG94 --clock --no-force --dir --gBGC --switch
