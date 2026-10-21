#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export ITDSEEKER_RSCRIPT=/tmp/flt3-itd-review-26cdfc1/dada2/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
sample=10808_hg38_RG
mkdir -p "review/runs/final_dada2/$sample"
python Nano_ITDseeker.py -b "bam_data/test_bam/$sample.bam" -o "review/runs/final_dada2/$sample" -s "$sample" -t 4 --haplotype-method dada2 --min-allele-frequency 0.01 > "review/runs/final_dada2/$sample/console.log" 2>&1
