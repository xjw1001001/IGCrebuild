#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YJR048W --paralog2 YEL039C --model HKY --no-clock --no-force --dir --gBGC
