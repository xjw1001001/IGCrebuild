#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YMR143W --paralog2 YDL083C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YMR143W --paralog2 YDL083C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YMR143W --paralog2 YDL083C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YMR143W --paralog2 YDL083C --model HKY --no-clock --no-force --dir --gBGC
