#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
export MAMBA_ROOT_PREFIX=/tmp/flt3-itd-review-26cdfc1/mamba
review/tools/bin/micromamba create --yes --prefix /tmp/flt3-itd-review-26cdfc1/dada2 --override-channels -c conda-forge -c bioconda bioconductor-dada2 bioconductor-shortread > review/dada2_native_setup.log 2>&1
