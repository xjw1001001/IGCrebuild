#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YBR191W --paralog2 YPL079W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR191W --paralog2 YPL079W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YBR191W --paralog2 YPL079W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR191W --paralog2 YPL079W --model HKY --no-clock --no-force --dir --gBGC
