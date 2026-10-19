#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
mkdir -p review/tools
curl --fail --location https://micro.mamba.pm/api/micromamba/linux-64/latest -o review/tools/micromamba.tar.bz2
 tar -xjf review/tools/micromamba.tar.bz2 -C review/tools bin/micromamba
export MAMBA_ROOT_PREFIX="$PWD/review/tools/mamba"
review/tools/bin/micromamba create --yes --prefix "$PWD/.validation-r" --override-channels -c conda-forge -c bioconda bioconductor-dada2 bioconductor-shortread > review/dada2_mamba_setup.log 2>&1
