#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YGL062W --paralog2 YBR218C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGL062W --paralog2 YBR218C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YGL062W --paralog2 YBR218C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGL062W --paralog2 YBR218C --model HKY --no-clock --no-force --dir --gBGC
