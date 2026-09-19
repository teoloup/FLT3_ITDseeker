# Reproducible environment for Nano_ITDseeker.
#
# Build:  docker build -t nano-itdseeker .
# Run:    docker run --rm -v "$PWD":/data nano-itdseeker \
#             -b /data/sample.bam -o /data/out -s SAMPLE -g hg38 -t 8 --html-report
#
# The image pins the external tools the pipeline shells out to (samtools) and
# the Python stack, including pymuscle5, which has no PyPI release and needs a
# C toolchain to build from source.

FROM python:3.10-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg

RUN apt-get update && apt-get install -y --no-install-recommends \
        samtools \
        build-essential \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/itdseeker

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip setuptools wheel \
    && pip install --no-cache-dir cython \
    && pip install --no-cache-dir -r requirements.txt

# --- haplotype-splitting backends -------------------------------------------
# Each is optional at runtime (--haplotype-method picks one), but the image
# carries all of them so a comparison run needs no extra setup.

# isONclust: the pip package is the Python implementation, which is enough here.
RUN pip install --no-cache-dir isONclust

# AmpliCI: source only, and it exits non-zero even on success, so the build
# checks that the binary exists rather than that it runs cleanly.
RUN git clone --depth 1 https://github.com/DormanLab/AmpliCI /opt/AmpliCI \
    && cd /opt/AmpliCI/src && cmake . && make -j"$(nproc)" \
    && ln -s /opt/AmpliCI/src/run_AmpliCI /usr/local/bin/run_AmpliCI \
    && test -x /usr/local/bin/run_AmpliCI

# dada2: R + Bioconductor. Bioconductor lags new R releases, so rather than
# pinning an R version here and having it drift, this uses a bioconda env, which
# ships dada2 with a compatible R alongside it. ITDSEEKER_RSCRIPT tells the
# dada2 backend which interpreter to use.
ENV MAMBA_ROOT_PREFIX=/opt/mamba \
    ITDSEEKER_RSCRIPT=/opt/dada2env/bin/Rscript
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C /usr/local bin/micromamba \
    && micromamba create -y -p /opt/dada2env -c conda-forge -c bioconda \
        bioconductor-dada2 bioconductor-shortread \
    && micromamba clean -y --all

COPY *.py dada2_cluster.R ./

# Fail the build rather than ship an image whose backends are silently missing.
RUN samtools --version | head -1 \
    && cutadapt --version \
    && python -c "import pymuscle5, pysam, Bio, sklearn; print('python deps OK')" \
    && isONclust --version \
    && test -x /usr/local/bin/run_AmpliCI \
    && "$ITDSEEKER_RSCRIPT" -e 'library(dada2); cat("dada2", as.character(packageVersion("dada2")), "OK\n")'

ENTRYPOINT ["python", "/opt/itdseeker/Nano_ITDseeker.py"]
