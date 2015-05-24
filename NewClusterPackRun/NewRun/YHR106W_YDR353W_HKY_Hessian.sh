#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YHR106W --paralog2 YDR353W --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YHR106W --paralog2 YDR353W --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YHR106W --paralog2 YDR353W --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YHR106W --paralog2 YDR353W --model HKY --no-clock --no-force --dir --gBGC
