#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export PATH=/tmp/flt3-medaka-review/venv/bin:/tmp/flt3-medaka-review/tools/bin:/usr/bin:/bin
medaka_consensus -h > review/runs/medaka/consensus_help.txt 2>&1
medaka tools list_models > review/runs/medaka/models.txt 2>&1
medaka tools resolve_model --auto_model consensus bam_data/test_bam/13697_2runs_hg38_RG.bam > review/runs/medaka/model_resolution.log 2>&1 || true
