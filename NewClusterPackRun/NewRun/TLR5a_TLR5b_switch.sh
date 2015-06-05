#!/bin/bash
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model HKY --clock --force --no-dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model HKY --clock --no-force --no-dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model HKY --clock --no-force --no-dir --gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model HKY --clock --no-force --dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model HKY --clock --no-force --dir --gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model MG94 --clock --force --no-dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model MG94 --clock --no-force --no-dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model MG94 --clock --no-force --no-dir --gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model MG94 --clock --no-force --dir --no-gBGC --switch
python Run_Fish.py --paralog1 TLR5a --paralog2 TLR5b --model MG94 --clock --no-force --dir --gBGC --switch
