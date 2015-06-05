#!/bin/bash
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model HKY --clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model HKY --clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model HKY --clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model HKY --clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model HKY --clock --no-force --dir --gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model MG94 --clock --force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model MG94 --clock --no-force --no-dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model MG94 --clock --no-force --no-dir --gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model MG94 --clock --no-force --dir --no-gBGC --switch
python Run_unfinished.py --paralog1 YPR159W --paralog2 YGR143W --model MG94 --clock --no-force --dir --gBGC --switch
