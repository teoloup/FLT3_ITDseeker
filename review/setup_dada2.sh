#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
/home/tloupis/miniconda3/bin/conda create --yes --prefix "$PWD/.validation-r" --override-channels -c conda-forge -c bioconda bioconductor-dada2 bioconductor-shortread > review/dada2_setup.log 2>&1
