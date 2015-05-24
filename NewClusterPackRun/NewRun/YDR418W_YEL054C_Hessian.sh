#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YDR418W --paralog2 YEL054C --model HKY --no-clock --no-force --dir --gBGC
