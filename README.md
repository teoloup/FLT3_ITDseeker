# Nano_ITDseeker (FLT3-ITD on Nanopore Reads)

`Nano_ITDseeker.py` detects and validates FLT3-ITD events from a BAM file of Nanopore reads.

## What The Program Does

1. Extract FLT3-region reads from BAM (`samtools`).
2. Trim amplicon with primer-aware `cutadapt` (`--rc` enabled).
3. Fit GMM on read lengths to identify WT/ITD peaks.
4. Optionally refine each detected ITD peak once into 2 local subpeaks.
5. Extract per-read insertion evidence by pairwise alignment.
6. Build per-peak consensus with MSA.
7. Build ITD synthetic references.
8. Validate reads by competitive alignment against WT + ITD refs.
9. Compute AF + strand-bias and export VCF (+ optional HTML report).

## Requirements

- Linux environment
- Python 3.9+ (recommended)
- External tools in `PATH`:
  - `samtools`
  - `cutadapt`
- Python packages:
  - `numpy`, `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, `scipy`
  - `biopython`, `pymuscle5`, `pysam`

## Input Requirements

- Coordinate-sorted BAM with index (`.bai`)
- Reads aligned to `hg38` or `hg19`
- FLT3 amplicon compatible with configured primers

## Quick Run

```bash
python Nano_ITDseeker.py \
  -b sample.bam \
  -o out \
  -s SAMPLE_01 \
  -g hg38 \
  -t 12 \
  --html-report
```

## Important CLI Options

- `-b, --bam`: input BAM (required)
- `-o, --output-folder`: output directory (required)
- `-s, --sample-name`: output prefix/sample id (required)
- `-g, --genome`: `hg38` or `hg19`
- `-t, --threads`: worker count
- `--min-allele-frequency`: minimum AF to report
- `--min-itd-size`, `--max-itd-size`: ITD size constraints
- `--per-peak-read-assignment-mode`: `manual|predict_proba|hybrid`
- `--force-number-of-peaks`: force initial GMM component count
- `--html-report`: write HTML report
- `--remove-intermediate-files`: delete `flt3_data` at end

Refinement options (new):
- `--disable-subpeak-refinement`
- `--min-reads-for-subpeak-refinement`
- `--min-subpeak-fraction`
- `--min-subpeak-distance`
- `--max-subpeak-sd`
- `--min-bic-gain-for-subpeak-split`

## Main Outputs

In output root:
- `<sample>_FLT3_ITD_calls.vcf`
- `<sample>_itd_report.html` (if `--html-report`)

In `flt3_data/`:
- `<sample>_itd_gmm_fit_plot.png`
- `<sample>_itd_insertions.tsv`
- `<sample>_itd_size_distribution.png`
- `<sample>_itd_consensus_seq.tsv`
- `<sample>_validation_refs.fasta`
- `<sample>_validation_read_support.tsv`
- per-ITD plots and alignments

## Notes On Strand Handling

- `cutadapt --rc` writes reverse-complemented reads when that orientation matches better.
- Reads are tagged with `strand` (`+` / `-`) for QC and downstream strand-bias stats.
- Alignment uses stored sequence directly (no dual-orientation re-search).

## Performance Tips

- **`-t 8` is the practical sweet spot.** The parallel stages scale well on
  their own (competitive validation measured 42.0s -> 4.7s, about 9x, going
  from 1 to 20 workers), but they are only 66-86% of single-thread runtime, so
  Amdahl caps total speedup at roughly 2.6-3.5x. Past 8 workers you buy about
  5% for 2.5x the cores.
- The remaining serial cost is dominated by the per-peak MUSCLE consensus and,
  behind it, the GMM component search.
- `--msa-max-unique` is the biggest single runtime knob, because MUSCLE is
  superlinear in panel size. The default of 150 was chosen against the
  validation set: raising it to 300 roughly triples the MSA stage (on one
  sample 9.1s -> 26.2s, total 35.6s -> 53.5s) without changing any call.
- Final validation is the heaviest parallel step; check logs for:
  - read count
  - batch count
  - elapsed time
  - comparisons/sec

## Install

See `requirements.txt`, or use the provided `Dockerfile`:

```bash
docker build -t nano-itdseeker .
docker run --rm -v "$PWD":/data nano-itdseeker \
    -b /data/sample.bam -o /data/out -s SAMPLE -g hg38 -t 8 --html-report
```

`pymuscle5` has no PyPI release and must be installed from source
(`pip install git+https://github.com/althonos/pymuscle5`), which needs Cython
and a C toolchain. `samtools` and `cutadapt` must both resolve on `PATH` --
installing cutadapt only inside a virtualenv that is not on `PATH` is not
enough, since the pipeline invokes it as a subprocess by name.

