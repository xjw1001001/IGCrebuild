#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YPL232W --paralog2 YMR183C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL232W --paralog2 YMR183C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YPL232W --paralog2 YMR183C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL232W --paralog2 YMR183C --model HKY --no-clock --no-force --dir --gBGC
