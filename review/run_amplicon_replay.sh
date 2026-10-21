#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH="$PWD/.validation-venv/bin:/usr/bin:/bin"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
mkdir -p review/runs/amplicon_sorter/pipeline_replay
python review/replay_amplicon_clusters.py > review/runs/amplicon_sorter/pipeline_replay/console.log 2>&1
