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

COPY *.py ./

# samtools and cutadapt must both resolve on PATH; fail the build if not.
RUN samtools --version | head -1 \
    && cutadapt --version \
    && python -c "import pymuscle5, pysam, Bio, sklearn; print('deps OK')"

ENTRYPOINT ["python", "/opt/itdseeker/Nano_ITDseeker.py"]
