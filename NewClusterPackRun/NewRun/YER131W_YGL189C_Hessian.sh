#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YER131W --paralog2 YGL189C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YER131W --paralog2 YGL189C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YER131W --paralog2 YGL189C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YER131W --paralog2 YGL189C --model HKY --no-clock --no-force --dir --gBGC
