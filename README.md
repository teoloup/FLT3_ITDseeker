# Nano_ITDseeker (FLT3-ITD on Nanopore Reads)

`Nano_ITDseeker.py` detects and validates FLT3-ITD events from a BAM file of Nanopore reads.

## What The Program Does

1. Extract FLT3-region reads from BAM (`samtools`).
2. Trim amplicon with primer-aware `cutadapt` (`--rc` enabled).
3. Fit GMM on read lengths to identify WT/ITD peaks.
4. Split each ITD peak into sequence haplotypes with DADA2 (default).
5. Extract per-read insertion evidence by pairwise alignment.
6. Build per-haplotype consensus with MSA, in WT reference context.
7. Build ITD synthetic references.
8. Validate reads by competitive alignment against WT + ITD refs.
9. Compute AF + strand-bias and export VCF (+ optional HTML report).

See `docs/ALGORITHM_LOGIC.md` for how each step works, and `docs/GUIDE.md`
(rendered as `docs/GUIDE.html`) for a user guide with worked examples.

## Repository Layout

```
Nano_ITDseeker.py   command-line entry point
itdseeker/          pipeline modules, plus dada2_cluster.R (the DADA2 wrapper)
docker/             Dockerfile; build with the repository root as context
docs/               ALGORITHM_LOGIC.md, GUIDE.md and the generated GUIDE.html
scripts/            simulate_itd_data.py, evaluate_haplotypes.py, build_guide_html.py
tests/              regression tests
benchmarks/         measurements behind the defaults
review/             validation records
```

## Haplotype Splitting With DADA2

Two different ITDs of the same size give reads of the same length, so the
length GMM puts them in one peak. Earlier versions tried to split peaks with a
second, local GMM on read length, but length cannot separate same-length ITDs.
That second GMM pass is no longer part of the pipeline. DADA2 now clusters each
peak's reads by sequence instead:

```
read lengths -> GMM peaks -> DADA2 per peak -> one candidate ITD per haplotype
             -> insertions, consensus, references, validation, VCF
```

- DADA2 runs on every ITD peak, not only peaks whose consensus looks mixed. An
  80/20 mixture yields a clean consensus with no `N`, so a mixed consensus is
  not a reliable trigger.
- An amplicon sequence variant (ASV) becomes its own haplotype only if it has at
  least `--min-haplotype-reads` reads and at least `--min-subpeak-fraction` of
  the peak. Reads from smaller ASVs are folded into the largest haplotype, so
  the AF denominator is unchanged.
- A split peak's haplotypes are named `<peak>_H1`, `<peak>_H2`, ... in
  decreasing read count.
- The WT peak is left alone unless `--cluster-wt-peak` is set.
- Separation is not guaranteed. A consensus that still contains unresolved bases
  is skipped with a warning rather than reported.

## Requirements

- Linux environment
- Python 3.9+ (recommended)
- External tools in `PATH`:
  - `samtools`
  - `cutadapt >= 5.2` (required for rightmost matching of the linked 3-prime primer)
- R with the Bioconductor packages `dada2` and `ShortRead`. Bioconductor often
  lags the newest R release, so if the system `Rscript` cannot load `dada2`, set
  `ITDSEEKER_RSCRIPT` to an Rscript that can, such as one from a bioconda env.
  The Docker image does this already.
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
- `--wt-peak-tolerance`: maximum WT peak offset from expected amplicon length (default 5 bp; must be smaller than `--min-itd-size`)
- `--html-report`: write HTML report
- `--remove-intermediate-files`: delete `flt3_data` at end

Haplotype options:
- `--haplotype-method dada2|none`: per-peak sequence clustering (default:
  `dada2`). `none` keeps the GMM length peaks unsplit.
- `--disable-subpeak-refinement`: same as `--haplotype-method none`, and
  overrides it.
- `--cluster-wt-peak`: also run DADA2 on the WT peak, to find ITDs hiding inside it
  (default: off)
- `--min-haplotype-reads`: minimum reads for an ASV to become a haplotype
  (default: 20). Peaks with fewer than twice this many reads are not clustered.
- `--min-subpeak-fraction`: minimum share of the parent peak for an ASV to
  become a haplotype (default: 0.15)
- `--dada2-omega-a`: DADA2 `OMEGA_A`, the p-value threshold for splitting off a
  new ASV. Lower is more conservative (default: `1e-40`).
- `--dada2-band-size`: DADA2 `BAND_SIZE` (default: 32, the value DADA2 documents
  for long, indel-prone reads)
- `--dada2-homopolymer-gap-penalty`: DADA2 `HOMOPOLYMER_GAP_PENALTY` (default:
  -1). Softens gaps in homopolymers, where ONT errors concentrate.

DADA2 is the only supported clustering tool. isONclust and AmpliCI have been
removed, and the GMM second pass is no longer used.

## Main Outputs

In output root:
- `<sample>_commands.jsonl` (exact external process arguments and stream routing)
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

In the temp directory (default `<output-folder>/temp`), under
`haplotypes/dada2_<peak>/`:
- `reads.fastq`, the peak's reads as given to DADA2
- `clusters.tsv`, each read's ASV
- `asvs.fasta`, the inferred ASV sequences with read counts

If the pipeline created the temp directory, it deletes it at the end of the run
unless `--log-level DEBUG` is set.

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

See `requirements.txt`, or use the provided `docker/Dockerfile`. Build from the
repository root, which is the build context:

```bash
docker build -f docker/Dockerfile -t nano-itdseeker .
docker run --rm -v "$PWD":/data nano-itdseeker \
    -b /data/sample.bam -o /data/out -s SAMPLE -g hg38 -t 8 --html-report
```

`pymuscle5` has no PyPI release and must be installed from source
(`pip install git+https://github.com/althonos/pymuscle5`), which needs Cython
and a C toolchain. `samtools` and `cutadapt` must both resolve on `PATH` --
installing cutadapt only inside a virtualenv that is not on `PATH` is not
enough, since the pipeline invokes it as a subprocess by name.


## Validation correctness

A WT length peak is optional. If none lies within `--wt-peak-tolerance` of
`--wt-amplicon-length`, candidate insertion sizes use the configured WT length.
When `--cluster-wt-peak` splits a WT peak, the closest eligible child keeps the
WT label; other children receive ITD aliases and undergo normal validation.

For multi-reference comparisons, a winning z-score must also have a softmax
probability margin of at least 0.05 over the runner-up. Exact and near ties remain
ambiguous. These softmax scores are relative alignment scores, not calibrated
variant probabilities. Breakpoint validation counts both insertions and deletions
within the ITD window against its gap allowance.

Insertion consensus is built with the WT reference context attached to each read's
insertion boundary. This keeps equivalent, shifted representations of a duplication
together instead of aligning rotated insertion payloads in isolation. The insertion
sequence and boundary are then recovered from the same consensus alignment.
`consensus_ins_pos_ref` is the resulting boundary; `raw_median_ins_pos_ref` preserves
the original read-position statistic. The legacy `median_ins_pos_ref` column now
aliases the consensus boundary for compatibility. MSA plots show the full contextual
allele; the consensus TSV and VCF continue to report the insertion payload.

The DADA2 backend honors `--threads` for both error learning and denoising.

Run the focused regression suite with the Python dependencies installed:

```bash
python -m unittest discover -s tests -v
```

## Exact external commands

Every samtools, Cutadapt, and DADA2 Rscript invocation is printed at INFO level
and appended to `<output-folder>/<sample>_commands.jsonl`. Records contain the
argument vector (`argv`), resolved executable, working directory, UTC timestamp,
POSIX-quoted display command, and stdin/stdout/stderr routing. Execution uses
`shell=False`; the argument vector is authoritative. The two samtools records
identify their pipe connection and FASTQ destination. Rscript arguments record
the wrapper, input/output paths, OMEGA_A, band size, homopolymer penalty, minimum
ASV reads, and threads, in that order; `itdseeker/dada2_cluster.R` contains the R calls.
Logs append across reruns and survive intermediate-file cleanup. They record
attempted invocations, not successful completion or tool-version provenance.
Replaying a DADA2 call requires its intermediate FASTQ to still exist.

Unresolved consensus bases are rejected before candidate validation. This avoids
reporting an ambiguous sequence as an ITD, but does not solve missing same-length
haplotypes or their effect on AF. See `review/CONSENSUS_VALIDATION.md`.
