#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YNL069C --paralog2 YIL133C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL069C --paralog2 YIL133C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YNL069C --paralog2 YIL133C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL069C --paralog2 YIL133C --model HKY --no-clock --no-force --dir --gBGC
