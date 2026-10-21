# Algorithm Logic And Concepts

This document explains the core workflow and the statistical/sequence-analysis principles used in `Nano_ITDseeker.py`.

## 0) Overview in plain language

An FLT3-ITD copies a stretch of the gene and inserts it right after the
original, so a read from an ITD allele is longer than a wild-type (WT) read by
exactly the size of the duplication. The pipeline works in two phases.

**Phase 1: propose candidates.** Several independent routes suggest which ITDs
might be present:

- the read-length distribution (one peak per ITD size, steps 2-4);
- sequence clustering inside each length peak, which separates different ITDs
  of the same or almost the same size (steps 3-4);
- a search for the same exact insertion recurring in many reads, which does not
  depend on the length peaks at all (step 5).

Each candidate gets a consensus sequence (step 6). Uncertain candidates are
dropped and duplicates are merged.

**Phase 2: test candidates.** Every candidate becomes a full reference sequence
(step 7). Every trimmed read is then aligned to WT and all candidates at once,
and goes to the reference it fits clearly best, if any (step 8). The AF of an
ITD is the share of reads it wins (step 9).

The phases are deliberately separate. Phase 1 may over-propose: a spurious
candidate wins few reads in phase 2 and is not reported. Phase 1 settings
(`--min-gmm-fraction`, `--min-haplotype-reads`, `--min-observed-haplotype-reads`)
decide what gets tested. Only `--min-allele-frequency` and `--min-itd-size`
decide what gets reported.

## 1) Pipeline Logic (Step by Step)

1. **Read extraction and primer trimming**
   - FLT3-region reads are extracted from BAM.
   - `cutadapt` trims linked primers with `--rc` enabled.
   - With `--rc`, cutadapt writes all output reads in one normalized orientation. Reads that were reverse-complemented by cutadapt are tagged with `rc` in the read header.
   - The pipeline stores:
     - `read_seq`: sequence exactly as output by cutadapt
     - `strand`: `+` or `-`, derived from the `rc` tag for QC and strand-bias statistics

2. **First-pass peak detection on read length**
   - The model fits read-length distribution with a GMM.
   - Components are filtered by quality constraints (minimum fraction, maximum SD).
   - Components that are too close can be merged to avoid over-fragmenting noise.
   - The component closest to the expected WT amplicon length is labeled `WT` if it lies within `--wt-peak-tolerance`; non-WT components become `ITD_*`.
   - A WT peak is optional. Without one, ITD sizes are measured against the configured `--wt-amplicon-length`.

3. **Per-peak haplotype splitting with DADA2**
   - Read length cannot separate two different ITDs of the same size, so they share one GMM peak. This step splits each peak by sequence.
   - It replaces the earlier second pass, which fitted a local `k=1` vs `k=2` length GMM inside each peak. That pass could only find close but different lengths, and it is no longer used.
   - Each ITD peak's primer-trimmed reads, with base qualities, are written to FASTQ and passed to `itdseeker/dada2_cluster.R`. The WT peak is included only with `--cluster-wt-peak`.
   - Peaks with fewer than `2 x --min-haplotype-reads` reads are not clustered.
   - DADA2 dereplicates the reads, learns an error model from that peak's own reads, and infers amplicon sequence variants (ASVs). Every read is then mapped back to its ASV.
   - Guardrails decide which ASVs become haplotypes:
     - An ASV needs at least `--min-haplotype-reads` reads (default 20) and at least `--min-subpeak-fraction` of the peak (default 0.15).
     - If fewer than two ASVs pass, the peak is left unsplit.
     - Reads from ASVs that fail, and reads DADA2 leaves unassigned, are folded into the largest haplotype rather than dropped. Dropping them would shrink the AF denominator and inflate every other call.
   - Haplotypes of peak `ITD_1` are named `ITD_1_H1`, `ITD_1_H2`, ... in decreasing read count. Each gets its own mean length, SD, read count and allele fraction, and is processed downstream as an independent candidate ITD.
   - When the WT peak is split, the child closest to the WT length (within `--wt-peak-tolerance`) keeps the `WT` label; the others become `ITD_WT_H*` candidates.
   - `--haplotype-method none` (or `--disable-subpeak-refinement`) skips this step and keeps the GMM peaks as they are.

4. **Per-peak insertion extraction, then a split by insertion length**
   - Reads assigned to each ITD peak or haplotype are aligned to the WT reference.
   - Insertions are inferred from alignment block structure.
   - Insertions are filtered to keep lengths compatible with the expected size range of that peak: its ITD size ±1.5 SD, but never narrower than ±`--min-gmm-peak-distance` (10 bp) while haplotype splitting is on. Peaks closer than that distance are merged, so a second ITD in the same peak can sit beyond 1.5 SD of its centre. Without the floor, ITDs 7–11 bp apart were lost here before the split below could see them. The floor stops halfway to the nearest other ITD peak, so neighbouring ranges never overlap.
   - DADA2's error model scores substitutions and tolerates indels, so two ITDs a few bp apart can stay in one peak. Their extracted insertion lengths separate cleanly, though, because only the inserted segment contributes errors.
   - Each peak is therefore split again when its insertion lengths form a second group with at least `--min-haplotype-reads` reads and `--min-subpeak-fraction` of the insertions, 3 bp or more away, and with the count between the groups dropping to at most half of the smaller one.
   - Reads go to the nearest length group, reads without an insertion to the largest, and the insertion table is relabelled with the new `_H*` aliases. Skipped with `--haplotype-method none`.

5. **Independent discovery of recurrent observed insertions**
   - This step runs on all trimmed reads before the length model, and does not use its peaks. It catches ITDs whose length peak is too small for the GMM (`--min-gmm-fraction`), ITDs hidden inside the WT peak, and sequences lost when a mixed peak gives an unresolved consensus.
   - Each read is located by two 20 bp anchor sequences taken from the WT amplicon, near its two ends. On the 336 bp reference they are at [15,35) and [305,325). An anchor may match with up to 3 substitutions but no indels. A read with no match, or with more than one possible anchor pair, is skipped.
   - The sequence between the anchors (the interior) is kept exactly as read. Anchor substitutions are replaced by WT and counted in the audit, so differences inside an anchor are not reported as variants.
   - Reads with an identical interior that is longer than WT are counted together. A group becomes a candidate only if it has:
     - at least `--min-observed-haplotype-reads` reads (default 20), with at least two on each strand;
     - median base quality of at least Q20;
     - a net size within `--min-itd-size` and `--max-itd-size`, and no unresolved bases;
     - an exact explanation as one pure insertion into WT: placing the inserted sequence into WT must reproduce the whole interior. Extra deletions, substitutions or scattered gaps disqualify the group, so sequencing errors never become short insertion calls.
   - How closely the insertion matches the adjacent WT copy is recorded as `DUP_IDENTITY` but is not a filter.
   - Accepted candidates are compared with the consensus candidates of step 6 and dropped when they give the same full allele. All candidates together are capped at `--max-itds-to-detect`.
   - Same-size candidates that differ by a single base at the junction are kept as separate sequences but flagged `SEQUENCE_UNCERTAINTY`: they may be one ITD written two ways, so their individual AFs are not clone sizes.
   - Every group and the reason it was accepted or rejected is written to `<sample>_observed_haplotype_audit.tsv`. `--disable-observed-rescue` turns this step off; the `--observed-anchor-*` options change anchor length, position and mismatch allowance.

6. **Per-peak MSA and consensus, in reference context**
   - Each read's insertion is placed back into the WT sequence at that read's insertion boundary, giving a full allele. Isolated insertion payloads can be cyclic rotations of the same duplication; aligning whole alleles keeps those together.
   - Alleles are deduplicated and weighted by abundance. The most abundant ones (up to `--msa-max-unique`, covering `--msa-min-weight-coverage` of reads) are aligned with MUSCLE5.
   - Only the varying part of each allele is aligned. In a tandem repeat an insertion can be written at a run of equivalent positions, and copies of one ITD share at least one of them even with sequencing errors in the insertion. The window is the narrowest span reaching such a position for 95% of the reads, plus 30 bp on each side, so a second ITD at another site widens it rather than being dropped. Outside it every allele is identical WT, so it is added back after alignment. This cuts MUSCLE's run time, which grows with sequence length. The few alleles that cannot place their insertion in the window are left out; at the real insertion they only carry gaps, which never vote on a base.
   - Weighted consensus is called per alignment column. Low-confidence columns are set to `N` using `--msa-base-threshold` and `--msa-min-col-coverage`.
   - A consensus containing any unresolved base is skipped with a warning, so an ambiguous sequence is never reported as an ITD.
   - Otherwise the consensus allele is realigned to WT. If it contains exactly one insertion within the size bounds, that insertion's sequence and position become the consensus payload and boundary (`consensus_ins_pos_ref`).
   - The largest per-column minority fraction is reported as `max_minor_fraction`. The pipeline warns when a consensus still looks mixed after clustering (`max_minor_fraction` > 0.15 or more than 5% `N`).
   - A per-ITD plot shows the insertion with 10 bp of WT each side: the per-base agreement, the consensus, and the most abundant aligned alleles with only their differences printed, so a recurrent second haplotype can be told apart from scattered sequencing errors.
   - Candidates whose consensus gives the same allele (insertion placed back into WT, so shifted representations match) are merged into the one with more reads. Identical references would make competitive validation call that ITD's reads ambiguous between them.

7. **ITD reference construction**
   - For each candidate, the consensus insertion is inserted into WT at the consensus boundary.
   - WT and all ITD references are written to a multi-reference FASTA for competitive validation.

8. **Competitive validation alignment**
   - Every primer-trimmed read is validated, not only those the GMM assigned to a peak. The GMM's job is to propose candidates; its narrow read windows (±2 SD, which a near-zero-width component can shrink to ±3 bp) dropped up to 12% of reads, more of them ITD than WT reads, and so biased AF low. `--validate-peak-reads-only` restores the earlier behaviour.
   - Each validation read is aligned to all references (WT + ITD references).
   - Scoring uses adjusted alignment score, softmax normalization, then best-vs-second-best comparison.
   - For multi-reference cases, z-score style normalization is used per read. A winner also needs a softmax margin of at least 0.05 over the runner-up; ties stay ambiguous.
   - Candidate ITD-supporting reads are additionally filtered by breakpoint/gap checks, which count both insertions and deletions in the ITD window.

9. **Final quantification and reporting**
   - ITD allele frequency is calculated from validated support reads.
   - Strand support is reported as per-ITD `plus,minus` counts.
   - Fisher exact test p-value for strand bias is reported in VCF INFO.
   - Only ITDs with AF of at least `--min-allele-frequency` and size of at least `--min-itd-size` are reported.
   - The main read-length plots are drawn at the end from the final read assignments, using each read's actual trimmed length: one panel with WT and all reported ITDs, and a zoom on the ITD reads. The earlier GMM fit is kept as a separate diagnostic.
   - For each reported ITD, up to 100 distinct supporting reads are aligned and shown as a consensus plot. This is for review only and never changes a call.
   - Output includes VCF and optional HTML report/plots.

## 2) Core Principles

### 2.1 Gaussian Mixture Model (GMM)

GMM assumes observed read lengths come from a mixture of latent Gaussian components:

\[
p(x) = \sum_{k=1}^{K} \pi_k \, \mathcal{N}(x \mid \mu_k, \sigma_k^2), \quad \sum_k \pi_k = 1
\]

- `mu_k` corresponds to expected read length mode.
- `sigma_k` reflects spread/noise around the mode.
- `pi_k` reflects component fraction (approximate abundance).

In this pipeline:
- WT reads form one dominant mode around WT amplicon length.
- ITDs form longer-length modes (often one per major ITD size family).
- Two ITDs of the same size form a single mode. The GMM cannot tell them apart; DADA2 (2.4) handles that.

### 2.2 EM algorithm used by GMM

GMM fitting uses Expectation-Maximization (EM):
- E-step: compute per-read responsibilities for each component.
- M-step: update component parameters (`pi`, `mu`, `sigma`) to maximize expected likelihood.
- Iterate until convergence.

Practical implication:
- EM can converge to local optima, so initialization and component constraints matter.

### 2.3 BIC for model selection

BIC balances goodness-of-fit against complexity:

\[
\mathrm{BIC} = k \ln(n) - 2\ln(\hat{L})
\]

- `k`: number of free parameters
- `n`: number of observations
- `L_hat`: maximized likelihood
- Lower BIC is preferred.

Pipeline use:
- Selects the number of components in the read-length GMM, unless `--force-number-of-peaks` is set.

### 2.4 DADA2 and amplicon sequence variants

DADA2 assumes every read is an error-containing copy of one of a small number of
true sequences (ASVs), and asks whether each distinct read sequence is too
abundant to be explained as errors from the ASVs already found.

- **Error model.** For an ASV `j` and a distinct read sequence `i`, `lambda_ji`
  is the probability that sequencing `j` produces `i`. It is the product of
  per-position transition probabilities (e.g. A->G at quality 12), learned from
  the data by `learnErrors`.
- **Abundance p-value.** If `j` has `n_j` reads, `i` is expected about
  `n_j * lambda_ji` times. With `a_i` observed copies, DADA2 computes the
  Poisson probability of seeing at least that many, given `i` was seen at all:

\[
p_A(j \rightarrow i) = \frac{1}{1 - \rho_{\mathrm{pois}}(n_j \lambda_{ji}, 0)} \sum_{a \ge a_i} \rho_{\mathrm{pois}}(n_j \lambda_{ji}, a)
\]

- **Splitting.** If that p-value falls below `OMEGA_A`, `i` becomes a new ASV and
  reads are re-partitioned. This repeats until no sequence is significant.
  Lower `OMEGA_A` makes a split harder; the default `1e-40` is deliberately
  conservative.

Why run it per peak:
- Within one peak every read comes from the same amplicon at nearly the same
  length, so one error model fits all of them. That is also the input DADA2 was
  built for: trimmed reads from a single amplicon.
- The error model is learned from that peak alone. If `learnErrors` fails, the
  wrapper falls back to DADA2's self-consistent estimate from the denoising step.

Adapting it to ONT reads:
- DADA2's error model covers substitutions. Indels only enter through the
  alignment between sequences, and ONT errors are mostly indels, especially in
  homopolymers.
- The wrapper therefore uses the settings DADA2 documents for long, indel-prone
  PacBio CCS reads: `BAND_SIZE = 32` (wider alignment band) and
  `HOMOPOLYMER_GAP_PENALTY = -1` (cheaper gaps inside homopolymers). Base
  qualities are used.

Why not trigger on `N`s:
- An `N` appears only when no base clears `--msa-base-threshold` (0.7), which
  needs the minority to exceed about 30%. At 80/20 the consensus is clean and the
  minor ITD is silently absorbed. DADA2 therefore runs on every ITD peak rather
  than in response to an ambiguous consensus.

Limits:
- Separation is not guaranteed. The mixed synthetic challenge with three 45 bp
  haplotypes in one peak is still not fully resolved; see
  `review/CONSENSUS_VALIDATION.md`.

### 2.5 Pairwise alignment for insertion detection

Each read is globally aligned to a target reference:
- Alignment blocks are parsed to detect where query advances more than target.
- These differences represent insertions in read relative to reference.
- We use score and percent identity, then positional checks near expected breakpoint for final validation.

### 2.6 MSA and weighted consensus

Per peak, alleles (the WT sequence carrying each read's insertion, trimmed to the window described in step 6) are aligned together:
- Weighted input emphasizes recurrent sequences and reduces influence of single-read noise.
- Consensus at each column is called by weighted support.
- Columns with low support/coverage are called as `N`.
- The insertion is recovered from the consensus allele afterwards, so its sequence and boundary come from the same alignment.

Why this matters:
- If two distinct ITDs are mixed into one peak, MSA columns become heterogeneous. A balanced mixture fills the consensus with `N`s and the candidate is skipped. An unbalanced one gives a clean majority consensus with a raised `max_minor_fraction`. DADA2 splitting (2.4) runs before this step to prevent both.

### 2.7 Softmax normalization and z-score in validation

For each read, alignment-to-reference scores are transformed with softmax:
- Produces a per-read relative support profile across references.
- Higher `beta` sharpens winner-take-most behavior.

Then:
- Binary case (WT vs one ITD): direct probability/delta interpretation.
- Multi-reference case: z-score-like normalization of per-read support values before best-hit calling, plus a minimum softmax margin of 0.05 between winner and runner-up.
- These scores are relative alignment support, not calibrated variant probabilities.

### 2.8 Strand bias testing

For each ITD:
- `SB` in FORMAT stores ITD-supporting counts as `plus,minus`.
- Fisher exact test compares ITD strand split against WT-supporting strand split.
- P-value is reported as `FISHER_P` in INFO.

## 3) Coordinate and VCF Representation Notes

- VCF `POS` is the genomic insertion anchor (1-based genomic coordinate as emitted by current pipeline mapping logic).
- VCF `REF` is anchored from WT sequence at the local insertion position.
- `ALT` stores the inserted ITD sequence payload.
- If `REF` falls back to `N`, this indicates missing/invalid local anchor metadata and should be treated as a data-path bug, not a biology signal.

## 4) Common Failure Modes and Why They Happen

- **Same-length ITDs in one peak**: the length GMM cannot separate them. DADA2 is meant to. If it does not, the consensus contains `N`s (candidate skipped) or shows a raised `max_minor_fraction` (minor ITD absorbed into the major one).
- **ITDs a few bp apart in one peak**: too close for the GMM, and DADA2 tolerates the length difference as indel noise. The insertion-length split (step 4) separates them when the minor ITD passes the guardrails. One below 20 reads or 15% of its peak is still absorbed into the major ITD, raising `max_minor_fraction`.
- **Over-splitting peaks**: DADA2 can call a sequencing-error variant a separate ASV. `OMEGA_A`, `--min-haplotype-reads` and `--min-subpeak-fraction` guard against this. A spurious haplotype usually collects little validated support and is not reported, but check it when it carries many reads.
- **Minor haplotype below the guardrails**: an ASV under 20 reads or 15% of its peak is folded into the largest haplotype, so a real low-level ITD sharing a length with a dominant one can be absorbed.
- **Low-level ITD missed by both routes**: an ITD below `--min-gmm-fraction` gets no length peak, and the recurrence search (step 5) needs 20 identical interiors with both strands and Q20. Noisy reads rarely repeat exactly, so a real ITD with few or error-rich reads can still be missed.
- **Poor validation separation**: references too similar or low-quality reads reduce best-vs-second-best margin.
- **Strand imbalance artifacts**: primer/trimming or read-quality effects can mimic biology; compare with WT split.

## 5) Main Runtime Data Objects

- `reads_df`
  - Core per-read table: `read_id`, `read_seq`, `read_qual`, `strand`, `read_len`, `gmm_peak_alias`, ambiguity flags. After splitting, `gmm_peak_alias` holds the haplotype alias.
- `comps`
  - Peak/component-level parameters: mean length, SD, fraction, aliases, AF-like metrics. Haplotypes carry `parent_peak_alias` and `is_refined_child`.
- `peak_subsets`
  - Mapping of peak or haplotype alias to supporting read IDs.
- `all_itd_insertions` / insertion tables
  - Per-read insertion evidence: insertion position, length, sequence, alignment stats.
- `df_cons`
  - Per-peak consensus summary: consensus insertion, `consensus_ins_pos_ref`, `raw_median_ins_pos_ref` and `max_minor_fraction`.
- `validation_results`
  - Competitive alignment metrics, support calls, filter reasons.
- `summary_df`
  - Final reported ITD-level AF and strand-bias statistics.

## 6) Useful References

- Cutadapt guide (`--revcomp` / `--rc` behavior):  
  https://cutadapt.readthedocs.io/en/latest/guide.html
- scikit-learn GMM user guide:  
  https://scikit-learn.org/stable/modules/mixture.html
- scikit-learn `GaussianMixture` API (`fit`, `bic`, `aic`):  
  https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html
- scikit-learn GMM model selection example (BIC/AIC):  
  https://sklearn.org/stable/auto_examples/mixture/plot_gmm_selection.html
- DADA2 documentation:  
  https://benjjneb.github.io/dada2/
- DADA2 paper (Callahan et al., Nature Methods, 2016):  
  https://doi.org/10.1038/nmeth.3869
- DADA2 on long PacBio CCS amplicons, source of the long-read settings (Callahan et al., Nucleic Acids Research, 2019):  
  https://doi.org/10.1093/nar/gkz569
- Biopython pairwise alignment tutorial (`Bio.Align.PairwiseAligner`):  
  https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html
- MUSCLE v5 documentation:  
  https://drive5.com/muscle5/manual/
- Original EM paper (Dempster, Laird, Rubin, 1977):  
  https://doi.org/10.1111/j.2517-6161.1977.tb01600.x
- Original BIC paper (Schwarz, 1978):  
  https://doi.org/10.1214/aos/1176344136
- Needleman-Wunsch global alignment paper (1970):  
  https://doi.org/10.1016/0022-2836(70)90057-4
