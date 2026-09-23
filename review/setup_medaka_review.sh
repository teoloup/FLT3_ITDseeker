#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/.."
mkdir -p review/runs/medaka
exec > review/runs/medaka/setup.log 2>&1
set -x
.validation-venv/bin/python -m venv /tmp/flt3-medaka-review/venv
/tmp/flt3-medaka-review/venv/bin/python -m pip install --upgrade pip
/tmp/flt3-medaka-review/venv/bin/python -m pip install 'torch==2.9.1+cpu' --index-url https://download.pytorch.org/whl/cpu
/tmp/flt3-medaka-review/venv/bin/python -m pip install 'medaka==2.2.2'
export MAMBA_ROOT_PREFIX=/tmp/flt3-medaka-review/mamba
review/tools/bin/micromamba create --yes --prefix /tmp/flt3-medaka-review/tools --override-channels -c conda-forge -c bioconda minimap2=2.30 samtools=1.22 htslib=1.22 bcftools=1.22
/tmp/flt3-medaka-review/venv/bin/medaka --version
