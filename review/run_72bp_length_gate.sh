#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
python review/benchmark_amplicon_sorter.py --sample 13697_2runs_hg38_RG --profile sc99 --native --length-diff-consensus 25
python review/benchmark_amplicon_sorter.py --sample sim_A --profile sc99 --native --length-diff-consensus 25
