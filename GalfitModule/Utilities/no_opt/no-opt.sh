#!bin/bash

python3 s2g_no-opt.py
cd sparcfire-out
for i in 123*; do galfit $i/"${i}.in"; done
cd ..
#sh model_png.sh
python3 panel_no-opt_models.py
