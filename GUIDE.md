# Nano_ITDseeker — user guide

Detecting FLT3 internal tandem duplications in Nanopore amplicon data: what the
tool does, how to run it, how to read what comes out, and where it can be wrong.

If you only want to run it, read [Quick start](#quick-start) and
[Reading the output](#reading-the-output). The rest explains the reasoning, with
worked examples from the validation set.

---

## Contents

1. [What an FLT3-ITD is, and why length is the starting point](#1-what-an-flt3-itd-is-and-why-length-is-the-starting-point)
2. [Quick start](#quick-start)
3. [How a sample flows through the pipeline](#3-how-a-sample-flows-through-the-pipeline)
4. [Worked example: a sample with three ITDs](#4-worked-example-a-sample-with-three-itds)
5. [The same-length problem](#5-the-same-length-problem)
6. [Reading the output](#reading-the-output)
7. [Parameters that actually matter](#7-parameters-that-actually-matter)
8. [Known limits](#8-known-limits)
9. [Validation](#9-validation)

---

## 1. What an FLT3-ITD is, and why length is the starting point

An FLT3-ITD duplicates a stretch of the juxtamembrane domain in tandem, so the
allele carries the duplicated bases inserted directly after the original copy.
Duplications in this cohort run from ~24 bp to ~190 bp.

The consequence the pipeline exploits: **an ITD read is longer than a wild-type
read by exactly the size of the duplication.** Amplify a fixed region, trim the
primers, and the read lengths fall into discrete peaks — one at the wild-type
amplicon length (336 bp by default), one per ITD, each offset by its size.

```
           WT (336bp)
              |
   reads      |            ITD +45bp            ITD +81bp
   per        ▇            (381bp)              (417bp)
   length    ▇▇▇              ▇                    ▇
            ▇▇▇▇▇            ▇▇▇                  ▇▇▇
   ─────────────────────────────────────────────────────── read length
```

Fitting a mixture model to those lengths gives you candidate ITDs, their sizes,
and a first estimate of how common each one is. That is the core of the method,
and for most samples it is sufficient.

Where it stops working is [section 5](#5-the-same-length-problem).

---

## Quick start

### Requirements

- Linux (or WSL). `samtools` and `cutadapt` must both resolve on `PATH` —
  installing cutadapt inside a virtualenv that is not on `PATH` is **not**
  enough, because the pipeline invokes it as a subprocess by name.
- Python 3.9+ and the packages in `requirements.txt`.
- `pymuscle5` has no PyPI release; install it from source:
  `pip install git+https://github.com/althonos/pymuscle5` (needs Cython and a C
  toolchain).

The `Dockerfile` does all of this, including the optional haplotype backends,
and fails the build if any of them is missing:

```bash
docker build -t nano-itdseeker .
docker run --rm -v "$PWD":/data nano-itdseeker \
    -b /data/sample.bam -o /data/out -s SAMPLE -g hg38 -t 8 --html-report
```

### Input

A coordinate-sorted, indexed BAM aligned to hg38 or hg19, containing the FLT3
amplicon. Reads may be in either orientation; the pipeline normalises them.

### A typical run

```bash
python Nano_ITDseeker.py \
    -b sample.bam \
    -o results/ \
    -s SAMPLE_01 \
    -g hg38 \
    -t 8 \
    --html-report
```

`-t 8` is deliberate — see [section 7](#7-parameters-that-actually-matter).

To report low-frequency ITDs, lower the reporting floor:

```bash
    --min-allele-frequency 0.01
```

The default of 0.05 will suppress real ITDs. In the validation set, two
fragment-analysis-confirmed ITDs sit at AF 0.040 and 0.017 — both are detected
and correctly quantified at the default, then filtered out before reporting.

---

## 3. How a sample flows through the pipeline

```
BAM
 │
 ├─ 1. extract region reads           samtools view + samtools fastq
 ├─ 2. trim primers                   cutadapt, linked adapters, --rc
 │                                    → reads normalised to one orientation,
 │                                      strand recorded per read
 ├─ 3. fit read lengths               Gaussian mixture; peaks = WT + candidates
 ├─ 4. split peaks into haplotypes    per-peak sequence clustering (DADA2)
 ├─ 5. extract insertions per peak    pairwise align each read to WT,
 │                                      read the inserted bases out
 ├─ 6. build a consensus per peak     weighted MSA over unique insertions
 ├─ 7. build one reference per ITD    WT with that consensus inserted
 ├─ 8. competitive validation         align every read against WT + all ITD
 │                                      references; each read votes once
 └─ 9. quantify and report            AF, strand bias, VCF, HTML
```

### Why step 8 exists

Steps 3–7 propose ITDs. Step 8 tests them. Every read is aligned against the
wild-type reference **and** every candidate ITD reference, and is assigned to
whichever it fits best. An ITD that was an artefact of read-length noise collects
no votes and disappears; a real one collects reads in proportion to its
abundance, which is where the reported allele frequency comes from.

This is also why AF is computed over *classified* reads rather than all reads:
a read that cannot be assigned confidently is not evidence for anything.

### Strand handling, briefly

`samtools fastq` restores each read to its original sequencing orientation, then
`cutadapt --rc` flips whichever reads need flipping and tags them. So the strand
recorded per read reflects the strand that was actually sequenced, which is what
makes the strand-bias test meaningful rather than circular.

---

## 4. Worked example: a sample with three ITDs

Sample 13697 from the validation set. Fragment analysis confirms three ITDs of
24 bp, 72 bp and 30 bp.

**Step 3 — read lengths.** The mixture fit finds a wild-type peak at ~336 bp and
two further peaks. The 24 bp and 30 bp ITDs are close enough in length that they
land in a single peak.

**Step 4 — splitting.** The length-based second pass tests that peak for
substructure and separates it into a 24 bp and a 30 bp component, because their
lengths genuinely differ.

**Steps 5–6 — insertions and consensus.** Each peak's reads are aligned to the
wild-type reference and the inserted bases extracted. For the 24 bp peak:

```
ITD_1_S1   len=24   n=531 reads   531 unique   0 N
           TCATATTCTCTGAAATCAACGTAG
```

Zero `N` is the signal that every read behind this peak agrees. If two different
ITDs had been mixed here, the disagreeing columns would come out as `N`.

**Step 8 — validation.** Four references now exist: WT, 24 bp, 72 bp, 30 bp.
Every read is aligned against all four.

**Step 9 — results.**

| ITD | size | position | AF | fragment analysis |
|---|---|---|---|---|
| ITD_1_S1 | 24 bp | chr13:28034125 | 8.7% | 24 bp, 8.5% |
| ITD_2 | 72 bp | chr13:28034103 | 8.0% | 72 bp, 10.5% |
| ITD_1_S2 | 30 bp | chr13:28034149 | 3.9% | 30 bp, 4.0% |

All three sizes match exactly. Note the third sits below the default
`--min-allele-frequency 0.05` and needs `0.01` to be reported.

---

## 5. The same-length problem

### The problem

Two different ITDs of the **same size** produce reads of the same length. They
are the same point in the feature the mixture model sees, so no amount of extra
length-based refinement can separate them. They end up in one peak, their reads
are pooled, and the consensus is built from a mixture of two sequences.

Where the two disagree, the consensus emits `N`:

```
one peak, three 45bp ITDs pooled:
  ANNNNTTTTCCAANNGNANNNNATNCTNCNGNAANN      ← 17 N in 36 bases
```

That output is useless, and worse, the reads supporting the second ITD are
counted towards the first, distorting both allele frequencies.

### What the pipeline does about it

The second pass clusters each peak's reads **by sequence**, not by length. The
flow is simply:

```
fit read lengths  ->  peaks
   for each peak  ->  subset its reads  ->  cluster them with DADA2
                                              -> one candidate ITD per cluster
   all candidates ->  insertions, consensus, references, validation, VCF
```

Two ITDs of the same length land in one length peak, and DADA2 separates them
there because it is looking at the sequences. Everything downstream is unchanged:
each cluster becomes a candidate ITD and is validated on the same footing as any
other.

On simulated data with three 45 bp ITDs sharing a peak this recovers two of them
plus the length-separable 30 bp one; the length-based pass recovered one. On all
five real validation samples, where the ITDs differ in size, it reproduces the
length-based result exactly.

**`--haplotype-method` selects the clustering step** if you want to compare:
`dada2` (default), `isonclust`, `amplici`, `gmm2pass` (the older length-based
split), `none`, or a chain like `gmm2pass+dada2`.

### Unbalanced mixtures leave no N

An `N` appears only when no base clears `--msa-base-threshold` (0.7), which
needs the minority to exceed about 30%. **A more lopsided mixture produces no `N`
at all** — at 80/20 the majority base clears the threshold and a minor ITD would
be silently absorbed. That is why the split happens unconditionally rather than
in response to `N`s.

The consensus table reports `max_minor_fraction`, the largest share of any column
disagreeing with the base called there, as a check that clustering worked:

| | max_minor_fraction |
|---|---|
| clean peaks, real samples (6 of them) | 0.024 – 0.101 |
| 80/20 mixture, before clustering | 0.217 |
| three-way same-length mixture, before clustering | 0.384 |

A value above ~0.15 in the output means a peak still looks mixed after
clustering, and those calls should be treated as provisional.

### Why dada2 and not the others

Three clustering backends are available. Measured against the validation set:

| method | 11531 (1 ITD) | 13697 (3 ITDs) | 14219 (2 ITDs) |
|---|---|---|---|
| `dada2` (default) | 1 ✓ | 3 ✓ | 2 ✓ |
| `gmm2pass` | 1 ✓ | 3 ✓ | 2 ✓ |
| `isonclust` | 1 ✓ | **2** ✗ | 2 ✓ |
| `amplici` | 1 ✓ | **2** ✗ | **1** ✗ |

`dada2` is the default because it is the only backend that matches the
length-based baseline everywhere while still separating same-length ITDs. `isonclust` merged the
30 bp ITD into the 24 bp peak; `amplici` fragmented real ITDs, in one case
splitting a 37%-AF ITD in two and then losing the call entirely. Neither is safe
to invoke automatically, though both remain available via `--haplotype-method`.

---

## Reading the output

```
results/
├── SAMPLE_FLT3_ITD_calls.vcf       ← the calls
├── SAMPLE_itd_report.html          ← the readable report
└── flt3_data/
    ├── SAMPLE_itd_gmm_fit_plot.png         read lengths and fitted peaks
    ├── SAMPLE_itd_size_distribution.png    insertion sizes across reads
    ├── SAMPLE_itd_insertions.tsv           per-read insertion evidence
    ├── SAMPLE_itd_consensus_seq.tsv        per-peak consensus (check N here)
    ├── SAMPLE_validation_refs.fasta        the references reads were tested against
    ├── SAMPLE_validation_read_support.tsv  per-read assignment and scores
    └── SAMPLE_<ITD>_*.png                  per-ITD plots
```

### The VCF

```
chr13  28034081  ITD_1  C  CCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATT  .  PASS
       TYPE=ITD;AF=0.41235;DP=13537;AF_GMM=0.40974;AF_FITTED=0.37672;
       FISHER_P=0.000478;ITD_LEN=45;INS_POS=28034081
       GT:DP:AF:AD:SB   0/1:13537:0.4124:7955,5582:2753,2829
```

| field | meaning |
|---|---|
| `POS` | 1-based plus-strand position of the base **before** the insertion |
| `REF` | that anchor base |
| `ALT` | the anchor base followed by the inserted sequence — standard VCF insertion form, so `bcftools norm` and VEP handle it correctly |
| `AF` | allele frequency from competitive validation. **This is the number to use.** |
| `AF_GMM` | allele frequency from read-length clustering, before validation |
| `AF_FITTED` | the mixture model's own weight for that peak |
| `DP` | classified reads contributing to AF |
| `AD` | reference-supporting, ITD-supporting |
| `SB` | ITD-supporting reads as `plus,minus` |
| `FISHER_P` | Fisher p for this ITD's strand split against the wild-type split. **Depth-sensitive — do not read alone.** |
| `STRAND_OR` | odds ratio for that comparison. 1.0 is no bias. This is the effect size. |
| `ITD_LEN` | duplication size in bp |

**Orientation.** FLT3 is transcribed from the minus strand of chr13. `REF` and
`ALT` are given on the **plus strand**, per VCF convention, so the sequence in
`ALT` is the reverse complement of the coding sequence.

**A caveat on `POS`.** An ITD sits in a tandem repeat, so a duplication of length
*L* inserted at position *P* is indistinguishable from the same duplication at
*P − L*; the aligner picks one. If you compare positions across tools or runs,
normalise first (`bcftools norm`) or compare modulo the ITD length.

### Three numbers worth checking

- **`AF` against `AF_GMM`.** They should be close. A large gap means competitive
  validation reassigned many reads, which is worth understanding before trusting
  the call.
- **`STRAND_OR` before `FISHER_P`.** The p-value scales with depth, so on a deep
  amplicon it flags differences far too small to matter. Sample 10808 in the
  validation set has an ITD at 49.3% plus against a wild type at 52.4% plus —
  three percentage points, odds ratio 0.885, no meaningful bias — and Fisher
  returns p = 0.0005 purely because there are 13,537 reads. The same split at 269
  reads gives p = 0.71. Read `STRAND_OR`: 1.0 is balanced, and it takes roughly a
  three-fold skew before a strand artefact is plausible. The report flags bias
  only when both a significant p and a three-fold odds skew are present, or when
  the variant is seen on essentially one strand.
- **`N` count in `*_itd_consensus_seq.tsv`.** Anything above zero means the reads
  behind that peak disagree. Clustering normally resolves this; if `N`s survive
  it, the peak may hold two ITDs the tool could not separate, and the call should
  be treated as provisional.

### The read-length plot

The plot shows every peak the model fitted, but a fitted peak is not a reported
ITD — competitive validation and the allele-frequency filter both sit downstream
and either can discard one. Peaks that did not survive are drawn dotted and grey
and labelled **not reported**, so the plot and the VCF agree.

Sample 11531 is the worked case: the length-based second pass split one peak into
components 3.18 bp apart, the smaller produced a 123 bp consensus against its own
expected 190.9 bp, and competitive validation gave it zero reads. It appears on
the plot, correctly marked, and not in the VCF. A discarded peak is usually a
spurious split rather than a missed ITD, but it is worth a look when it carries
many reads.

### The HTML report

One card per ITD, carrying the size, genomic position, allele frequency with its
read counts, strand balance with the Fisher test, and the inserted sequence
itself (ambiguous bases highlighted). Whether the frame is preserved is shown,
since `ITD_LEN % 3 != 0` means a frameshift rather than an in-frame duplication.
Everything is inlined, so the file can be archived or emailed on its own.

---

## 7. Parameters that actually matter

Defaults were chosen against the validation set and the simulator; the
measurements are in [`benchmarks/`](benchmarks/).

| option | default | why |
|---|---|---|
| `--min-allele-frequency` | 0.05 | **Often too high.** Real fragment-analysis-confirmed ITDs in this cohort sit at 0.040 and 0.017. Use `0.01` unless you need the stricter floor. |
| `-t, --threads` | 1 | **8 is the practical maximum.** The parallel stages scale ~9×, but they are only 66–86% of runtime, so total speedup caps near 2.6–3.5×. Past 8 workers you gain ~5% for 2.5× the cores. |
| `--haplotype-method` | `dada2` | Per-peak sequence clustering. Separates same-length ITDs; matches the length-based result on every real sample. |
| `--msa-max-unique` | 150 | Largest single runtime lever — MUSCLE is superlinear in panel size. Raising to 300 roughly triples the consensus stage without changing any call. |
| `--min-itd-size` / `--max-itd-size` | 12 / 300 | Hard bounds on reportable duplication size. |
| `--wt-amplicon-length` | 336 | Must match your amplicon. It defines which peak is wild-type, and therefore every ITD size. |

If you change the amplicon or primers, `--wt-amplicon-length` and the primer
sequences in `bam_extractor.py` must both be updated, and the reference amplicon
in `Nano_ITDseeker.py` with them.

---

## 8. Known limits

- **ITDs identical in both length and position**, differing only by substituted
  bases, are not separated. Only `amplici` resolved this in simulation, and it
  is too aggressive on real data to use by default.
  This pattern is biologically uncommon — co-occurring ITDs usually differ in
  size or breakpoint — but a surviving `N` in the consensus is the sign of it.
- **Detection depends on read-length clustering.** An ITD whose reads do not form
  a distinct length peak, or that falls below `--min-gmm-fraction` (default 1%),
  will not be proposed and therefore cannot be validated. Reads assigned to the
  wild-type peak are never examined for insertions.
- **Very small ITDs** close to the wild-type length can be absorbed into the
  wild-type peak.
- **`POS` is not left-normalised.** See the caveat above.
- **Genotype is always reported `0/1`.** These are somatic calls; the field is a
  placeholder, not a zygosity determination.
- **The `AD` reference count** is depth minus this ITD's reads, so in a
  multi-ITD sample it includes reads supporting the other ITDs.

---

## 9. Validation

Five samples, seven fragment-analysis-confirmed ITDs, one ITD-negative control.

| sample | ITDs | sizes | result |
|---|---|---|---|
| 10808 | 1 | 45 bp | ✓ size and position exact |
| 11531 | 1 | 189 bp | ✓ |
| 13697 | 3 | 24, 72, 30 bp | ✓ all three |
| 14219 | 2 | 81, 30 bp | ✓ both |
| 14417 | 0 | — | ✓ no calls |

**7/7 detected, 0 false positives.** Every ITD length and genomic position
matches the previous validated version of the tool exactly, and the inserted
sequences are byte-identical. Allele frequencies agree with fragment analysis to
within the difference expected between the two assays.

Measured read accuracy on this cohort, from the ITD-negative control: median
99.4% identity, with indels running about 1.3× substitutions — consistent with
R10.4 chemistry and a recent super-accuracy basecalling model.

To regenerate or extend the evidence:

```bash
# build a synthetic sample with known haplotypes
python simulate_itd_data.py --scenario A --reads 8000 --out-dir sim_data

# run it
python Nano_ITDseeker.py -b sim_data/sim_A.bam -o sim_out -s simA \
    -g hg38 -t 8 --min-allele-frequency 0.01

# score the calls against the truth tables
python evaluate_haplotypes.py --run-dir sim_out --sample simA \
    --truth-dir sim_data --scenario A --method gmm2pass
```

`simulate_itd_data.py` models the error profile measured from the negative
control, and scenario A deliberately includes ITDs that share a length, so it
exercises the failure mode that motivates per-peak sequence clustering. `evaluate_haplotypes.py`
reports recovery, over-splitting, allele-frequency error, ambiguous-base counts
and read-level clustering agreement.
