#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YGL133W --paralog2 YPL216W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGL133W --paralog2 YPL216W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YGL133W --paralog2 YPL216W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGL133W --paralog2 YPL216W --model HKY --no-clock --no-force --dir --gBGC
