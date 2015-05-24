#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YJL177W --paralog2 YKL180W --model MG94 --no-clock --no-force --dir --gBGC
