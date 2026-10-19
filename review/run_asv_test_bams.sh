#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export ITDSEEKER_RSCRIPT=/tmp/flt3-itd-review-26cdfc1/dada2/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=4
"$ITDSEEKER_RSCRIPT" -e 'library(dada2); library(ShortRead); cat("DADA2",as.character(packageVersion("dada2")),"\n")'
python review/run_test_bams.py dada2
