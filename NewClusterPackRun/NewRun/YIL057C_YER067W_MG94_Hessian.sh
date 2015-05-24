#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YIL057C --paralog2 YER067W --model MG94 --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YIL057C --paralog2 YER067W --model MG94 --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YIL057C --paralog2 YER067W --model MG94 --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YIL057C --paralog2 YER067W --model MG94 --no-clock --no-force --dir --gBGC
