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

- Use `-t` close to available physical cores.
- Final validation is the heaviest step; check logs for:
  - read count
  - batch count
  - elapsed time
  - comparisons/sec

