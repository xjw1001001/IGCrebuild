#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YAL056W --paralog2 YOR371C --model HKY --no-clock --no-force --dir --gBGC
