#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YGR043C --paralog2 YLR354C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGR043C --paralog2 YLR354C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YGR043C --paralog2 YLR354C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YGR043C --paralog2 YLR354C --model HKY --no-clock --no-force --dir --gBGC
