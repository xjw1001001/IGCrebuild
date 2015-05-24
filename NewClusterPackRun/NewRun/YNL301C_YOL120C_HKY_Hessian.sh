#!/bin/bash
python GenerateInfomationMatrix.py --paralog1 YNL301C --paralog2 YOL120C --model HKY --no-clock --no-force --no-dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL301C --paralog2 YOL120C --model HKY --no-clock --no-force --no-dir --gBGC
python GenerateInfomationMatrix.py --paralog1 YNL301C --paralog2 YOL120C --model HKY --no-clock --no-force --dir --no-gBGC
python GenerateInfomationMatrix.py --paralog1 YNL301C --paralog2 YOL120C --model HKY --no-clock --no-force --dir --gBGC
