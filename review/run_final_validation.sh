#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export ITDSEEKER_RSCRIPT=/tmp/flt3-itd-review-26cdfc1/dada2/bin/Rscript
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m unittest discover -s tests -v > review/runs/final_unit_tests.log 2>&1
python Nano_ITDseeker.py --help > review/runs/final_help.txt
for sample in sim_A sim_no_wt; do
    mkdir -p "review/runs/final_dada2/$sample"
    python Nano_ITDseeker.py -b "review/synthetic/$sample.bam" -o "review/runs/final_dada2/$sample" -s "$sample" -t 4 --haplotype-method dada2 --min-allele-frequency 0.01 > "review/runs/final_dada2/$sample/console.log" 2>&1
done
