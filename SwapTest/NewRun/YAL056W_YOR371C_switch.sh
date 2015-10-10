#!/bin/bash
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model MG94 --no-clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model MG94 --no-clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model MG94 --no-clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model MG94 --no-clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YAL056W --paralog2 YOR371C --model MG94 --no-clock --no-force --dir --gBGC --switch
