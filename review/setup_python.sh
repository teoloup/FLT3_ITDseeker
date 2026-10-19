#!/usr/bin/env bash
set -eu
cd /mnt/c/Users/HemaLab/Desktop/CODING/FLT3_TOOL/FLT3_ITDseeker
/home/tloupis/miniconda3/envs/cutadapt/bin/python -m venv --system-site-packages .validation-venv
.validation-venv/bin/python -m pip install --upgrade cutadapt
.validation-venv/bin/cutadapt --version
