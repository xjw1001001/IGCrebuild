#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model HKY --no-clock --no-force --dir --gBGC
