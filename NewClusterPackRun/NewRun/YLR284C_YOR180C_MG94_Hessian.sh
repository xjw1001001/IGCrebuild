#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YLR284C --paralog2 YOR180C --model MG94 --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YLR284C --paralog2 YOR180C --model MG94 --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YLR284C --paralog2 YOR180C --model MG94 --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YLR284C --paralog2 YOR180C --model MG94 --no-clock --no-force --dir --gBGC
