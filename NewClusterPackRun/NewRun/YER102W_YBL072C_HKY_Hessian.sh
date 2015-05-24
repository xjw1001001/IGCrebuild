#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YER102W --paralog2 YBL072C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YER102W --paralog2 YBL072C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YER102W --paralog2 YBL072C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YER102W --paralog2 YBL072C --model HKY --no-clock --no-force --dir --gBGC
