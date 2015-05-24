#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YPL037C --paralog2 YDR252W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL037C --paralog2 YDR252W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YPL037C --paralog2 YDR252W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL037C --paralog2 YDR252W --model HKY --no-clock --no-force --dir --gBGC
