#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YIR033W --paralog2 YKL020C --model MG94 --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YIR033W --paralog2 YKL020C --model MG94 --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YIR033W --paralog2 YKL020C --model MG94 --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YIR033W --paralog2 YKL020C --model MG94 --no-clock --no-force --dir --gBGC
