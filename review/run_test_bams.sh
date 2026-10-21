#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH=/mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker/.validation-venv/bin:/usr/bin:/bin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python review/run_test_bams.py prepare
python review/run_test_bams.py gmm2pass
