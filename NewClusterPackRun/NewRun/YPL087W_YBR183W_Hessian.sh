#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YPL087W --paralog2 YBR183W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL087W --paralog2 YBR183W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YPL087W --paralog2 YBR183W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YPL087W --paralog2 YBR183W --model HKY --no-clock --no-force --dir --gBGC
