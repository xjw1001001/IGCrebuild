#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YDR518W --paralog2 YCL043C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YDR518W --paralog2 YCL043C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YDR518W --paralog2 YCL043C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YDR518W --paralog2 YCL043C --model HKY --no-clock --no-force --dir --gBGC
