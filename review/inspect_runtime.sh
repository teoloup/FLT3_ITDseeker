#!/usr/bin/env bash
set -eu
/home/tloupis/miniconda3/bin/conda env list
for env in /home/tloupis/miniconda3 /home/tloupis/miniconda3/envs/pyRjup /home/tloupis/miniconda3/envs/cutadapt; do
  echo "$env"
  "$env/bin/python" -c 'import importlib.util; print({x: bool(importlib.util.find_spec(x)) for x in ["numpy","pandas","scipy","sklearn","Bio","pysam","pymuscle5","cutadapt","matplotlib","seaborn"]})' || true
done
Rscript -e 'cat("dada2 available: ",requireNamespace("dada2",quietly=TRUE),"\n")'
