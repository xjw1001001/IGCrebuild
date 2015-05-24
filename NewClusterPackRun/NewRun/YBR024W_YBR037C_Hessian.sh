#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YBR024W --paralog2 YBR037C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR024W --paralog2 YBR037C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YBR024W --paralog2 YBR037C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YBR024W --paralog2 YBR037C --model HKY --no-clock --no-force --dir --gBGC
