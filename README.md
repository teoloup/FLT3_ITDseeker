# Nano_ITDseeker

Nano_ITDseeker finds FLT3 internal tandem duplications (FLT3-ITDs) in Nanopore
amplicon sequencing data. It takes an aligned BAM file, works out which ITDs are
present, checks each one read by read, and reports its size, sequence, position
and allele frequency (AF) in a VCF file and, if requested, an HTML report.

## How it works, in plain language

An FLT3-ITD is a piece of the gene copied and inserted right after the original.
A read from an ITD allele is therefore longer than a wild-type (WT) read by
exactly the size of the duplication. After the primers are trimmed, WT reads all
have about the same length (336 bp), and each ITD shows up as a separate group of
longer reads.

The pipeline uses that in two phases: first it **proposes** candidate ITDs, then
it **tests** every candidate against all the reads.

**Proposing candidates**

1. **Collect the reads.** Reads covering the FLT3 amplicon are taken from the
   BAM, and the primers are trimmed off. All reads are put in the same
   orientation, and the strand each came from is kept.
2. **Group reads by length.** A statistical model (a Gaussian mixture model)
   finds the length peaks: one for WT, and one per ITD size.
3. **Separate ITDs that share a length.** Two different ITDs of the same size
   fall in the same length peak. Each peak is sorted by sequence with DADA2, a
   tool that tells real sequence variants apart from sequencing errors. Peaks
   holding two ITDs only a few bases apart are also split, using the size of the
   insertion found in each read.
4. **Look for ITDs the length model missed.** Independently of the length peaks,
   the pipeline looks for the exact same insertion repeated in many reads (by
   default at least 20, from both strands and of good quality). This catches
   small ITDs whose length peak is too weak to be picked up.
5. **Build each candidate's sequence.** For each candidate, the reads are
   aligned together and a consensus is called base by base. Candidates whose
   sequence is still uncertain are dropped rather than reported, and candidates
   that turn out to be the same ITD are merged.

**Testing candidates**

6. **Make reference sequences.** Each candidate's insertion is placed into the
   WT sequence to make a full ITD reference.
7. **Let the references compete for each read.** Every trimmed read is aligned
   to WT and to all ITD references. A read counts for an ITD only if it matches
   that reference clearly better than all others and the ITD junction itself
   looks correct. Reads that fit two references equally well count for neither.
8. **Report.** AF is the share of reads won by each ITD. Strand balance is
   checked against WT reads. ITDs at or above the minimum AF (5% by default) are
   written to the VCF and report.

Proposing and testing are kept separate on purpose. A spurious candidate only
costs some extra alignment time, because it wins few reads and is not reported.
The final AF always comes from the read-by-read test, not from the length peaks.

For the full method, see [docs/ALGORITHM_LOGIC.md](docs/ALGORITHM_LOGIC.md). A
longer user guide is in [docs/GUIDE.md](docs/GUIDE.md) (also as
[docs/GUIDE.html](docs/GUIDE.html)).

## Input

- A coordinate-sorted, indexed BAM file (`.bam` with `.bai`).
- Reads aligned to `hg38` or `hg19`.
- An FLT3 amplicon made with the primers the pipeline expects (WT length 336 bp
  after trimming).

## Installation

The pipeline runs on Linux. On Windows, use WSL or Docker.

### Option 1: Docker (recommended)

The Docker image contains every dependency, including R and DADA2. Build it
from the repository root:

```bash
git clone https://github.com/teoloup/FLT3_ITDseeker.git
cd FLT3_ITDseeker
docker build -f docker/Dockerfile -t nano-itdseeker .
```

The image can also be converted for Apptainer/Singularity, where it provides a
`Nano_ITDseeker` command.

### Option 2: Manual installation

You need:

- Python 3.9 or newer
- `samtools`
- `cutadapt` 5.2 or newer
- R with the Bioconductor packages `dada2` and `ShortRead`
- a C compiler (to build `pymuscle5`)

One way to get the tools is with conda/mamba:

```bash
mamba create -n itdseeker -c conda-forge -c bioconda \
    python=3.10 samtools "cutadapt>=5.2" \
    bioconductor-dada2 bioconductor-shortread compilers
mamba activate itdseeker

git clone https://github.com/teoloup/FLT3_ITDseeker.git
cd FLT3_ITDseeker
pip install cython
pip install -r requirements.txt
```

Notes:

- `pymuscle5` has no PyPI release. `requirements.txt` installs it from GitHub,
  which needs Cython and a C compiler.
- `samtools` and `cutadapt` must be on `PATH`, because the pipeline calls them
  by name.
- If the default `Rscript` cannot load `dada2`, point the pipeline at one that
  can: `export ITDSEEKER_RSCRIPT=/path/to/env/bin/Rscript`.

Check the installation by running the tests (`MPLBACKEND=Agg` lets the plot
tests run on a machine without a display):

```bash
MPLBACKEND=Agg python -m unittest discover -s tests -v
```

## Running it

Minimal run:

```bash
python Nano_ITDseeker.py -b sample.bam -o results -s SAMPLE -g hg38 -t 8 --html-report
```

With Docker, mount the folder that holds your data:

```bash
docker run --rm -v "$PWD":/data nano-itdseeker \
    -b /data/sample.bam -o /data/results -s SAMPLE -g hg38 -t 8 --html-report
```

`-t 8` is a good choice. More threads help little, because part of the pipeline
runs on a single thread.

Run `python Nano_ITDseeker.py --help` to see every option.

### Main options

| Option | Meaning | Default |
|---|---|---|
| `-b, --bam` | Input BAM file | required |
| `-o, --output-folder` | Output folder | required |
| `-s, --sample-name` | Sample name, used as the prefix of every output file | required |
| `-g, --genome` | Reference genome of the BAM: `hg38` or `hg19` | `hg38` |
| `-t, --threads` | Number of threads | 1 |
| `--html-report` | Also write an HTML report | off |
| `--min-allele-frequency` | Lowest AF to report | 0.05 |
| `--min-itd-size`, `--max-itd-size` | Allowed ITD sizes (bp) | 12, 300 |
| `--max-itds-to-detect` | Maximum number of candidate ITDs tested | 6 |
| `-q, --min-mapping-quality` | Minimum mapping quality of extracted reads | 20 |
| `--remove-intermediate-files` | Delete the `flt3_data/` folder at the end | off |
| `--log-level` | `INFO` or `DEBUG` (`DEBUG` also keeps the temp folder) | `INFO` |

### Options for finding candidates

These change which candidates are proposed. They do not change how a
candidate's AF is measured.

| Option | Meaning | Default |
|---|---|---|
| `--min-gmm-fraction` | Smallest share of reads for a length peak to be kept | 0.01 |
| `--min-gmm-peak-distance` | Length peaks closer than this (bp) are merged | 10 |
| `--force-number-of-peaks` | Fix the number of length peaks instead of choosing it automatically | automatic |
| `--haplotype-method` | `dada2` to split peaks by sequence, `none` to keep them as they are | `dada2` |
| `--min-haplotype-reads` | Reads needed for a sequence group to become its own candidate | 20 |
| `--min-subpeak-fraction` | Share of its peak needed for a sequence group to become its own candidate | 0.15 |
| `--cluster-wt-peak` | Also split the WT peak by sequence, to find ITDs hidden inside it | off |
| `--disable-observed-rescue` | Turn off the search for repeated exact insertions (step 4) | on |
| `--min-observed-haplotype-reads` | Identical reads needed to propose a candidate in step 4 | 20 |
| `--validate-peak-reads-only` | Test only reads inside a length peak, instead of every read | off |

`--min-gmm-fraction` and `--min-allele-frequency` are independent on purpose:
the first decides which peaks become candidates, the second which tested ITDs
are reported.

Advanced settings for DADA2 (`--dada2-*`), the consensus (`--msa-*`) and the
step-4 search (`--observed-anchor-*`) are listed by `--help` and explained in
[docs/ALGORITHM_LOGIC.md](docs/ALGORITHM_LOGIC.md).

## Output

In the output folder:

| File | Content |
|---|---|
| `<sample>_FLT3_ITD_calls.vcf` | Reported ITDs: position, inserted sequence, size, AF, supporting reads per strand, strand-bias test |
| `<sample>_itd_report.html` | Report with read-length plots, each ITD's sequence and support (with `--html-report`) |
| `<sample>_commands.jsonl` | Exact `samtools`, `cutadapt` and DADA2 commands that were run |

In `flt3_data/` (intermediate files and plots):

| File | Content |
|---|---|
| `<sample>_itd_consensus_seq.tsv` | Every candidate and its consensus sequence |
| `<sample>_itd_insertions.tsv` | The insertion found in each read |
| `<sample>_validation_refs.fasta` | WT and ITD references used for testing |
| `<sample>_validation_read_support.tsv` | Which reference each read supported, and why |
| `<sample>_observed_haplotype_audit.tsv` | Candidates from step 4 and why each was accepted or rejected |
| `<sample>_itd_gmm_fit_plot.png` | Read lengths of WT and reported ITD reads |
| `<sample>_<ITD>_final_support_MSA_consensus.png` | Alignment of each reported ITD's supporting reads |

A temp folder (default `<output-folder>/temp`) holds the DADA2 inputs and
results. It is deleted at the end unless `--log-level DEBUG` is set.

### Reading the results

- **AF** is the share of tested reads supporting the ITD.
- **Strand bias**: `SB` gives the ITD's plus/minus read counts. `FISHER_P` tests
  whether that split differs from WT's, and `STRAND_OR` says by how much (1.0
  means no bias). With deep coverage a small, harmless difference can still give
  a low p-value, so look at both. A strong bias can point to an artefact.
- **`SEQUENCE_UNCERTAINTY`** marks ITDs of the same size whose sequences differ
  by one base at the junction. They may be one ITD described two ways rather
  than two clones, so their separate AFs should not be read as clone sizes.

## Limits

- Separating ITDs is not guaranteed. A minor ITD with the same length as a major
  one can be absorbed into it when it has fewer than 20 reads or less than 15%
  of its peak.
- Candidates with an uncertain consensus are dropped rather than reported, so a
  real but poorly resolved ITD can be missed.
- Only ITDs at or above `--min-allele-frequency` and `--min-itd-size` are
  reported.

## Repository layout

```
Nano_ITDseeker.py   command-line entry point
itdseeker/          pipeline modules, plus dada2_cluster.R (the DADA2 wrapper)
docker/             Dockerfile (build with the repository root as context)
docs/               method description and user guide
scripts/            simulation, evaluation and guide-building helpers
tests/              regression tests
benchmarks/         measurements behind the defaults
review/             validation records
```
