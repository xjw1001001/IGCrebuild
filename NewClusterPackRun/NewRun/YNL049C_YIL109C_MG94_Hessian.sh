#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model MG94 --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model MG94 --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model MG94 --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL049C --paralog2 YIL109C --model MG94 --no-clock --no-force --dir --gBGC
