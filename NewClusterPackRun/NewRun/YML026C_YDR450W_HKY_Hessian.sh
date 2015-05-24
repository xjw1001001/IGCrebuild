#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YML026C --paralog2 YDR450W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YML026C --paralog2 YDR450W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YML026C --paralog2 YDR450W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YML026C --paralog2 YDR450W --model HKY --no-clock --no-force --dir --gBGC
