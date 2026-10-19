#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export ITDSEEKER_RSCRIPT=/tmp/flt3-itd-review-26cdfc1/dada2/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python review/run_clean_control.py
