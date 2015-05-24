#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YBR117C --paralog2 YPR074C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR117C --paralog2 YPR074C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YBR117C --paralog2 YPR074C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR117C --paralog2 YPR074C --model HKY --no-clock --no-force --dir --gBGC
