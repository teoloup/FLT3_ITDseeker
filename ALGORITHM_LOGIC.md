# Algorithm Logic And Concepts

This document explains the core workflow and the statistical/sequence-analysis principles used in `Nano_ITDseeker.py`.

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
   - The component closest to expected WT amplicon length is labeled `WT`; non-WT components become `ITD_*`.

3. **One-level local subpeak refinement**
   - For each first-pass ITD peak, the pipeline tests whether that peak should be split into two local subpeaks.
   - This uses a local `k=1` vs `k=2` comparison with BIC and additional guardrails (minimum child size, separation, SD limits).
   - Refinement is intentionally limited to one level to stay conservative and avoid overfitting.

4. **Per-peak insertion extraction**
   - Reads assigned to each ITD peak are aligned to WT reference.
   - Insertions are inferred from alignment block structure.
   - Insertions are filtered to keep lengths compatible with the expected size range of that peak.

5. **Per-peak MSA and consensus**
   - Insertion sequences are deduplicated and weighted by abundance.
   - MSA is computed per peak.
   - Weighted consensus is called per alignment column.
   - Low-confidence columns are set to `N` using thresholds.

6. **ITD reference construction**
   - For each peak, consensus insertion is inserted into WT at the median insertion position.
   - WT and all ITD references are written to a multi-reference FASTA for competitive validation.

7. **Competitive validation alignment**
   - Each validation read is aligned to all references (WT + ITD references).
   - Scoring uses adjusted alignment score, softmax normalization, then best-vs-second-best comparison.
   - For multi-reference cases, z-score style normalization is used per read.
   - Candidate ITD-supporting reads are additionally filtered by breakpoint/gap checks.

8. **Final quantification and reporting**
   - ITD allele frequency is calculated from validated support reads.
   - Strand support is reported as per-ITD `plus,minus` counts.
   - Fisher exact test p-value for strand bias is reported in VCF INFO.
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
- First pass: supports component-count selection.
- Local refinement: split accepted only when `BIC(k=1) - BIC(k=2)` exceeds threshold and biologic/quality constraints also pass.

### 2.4 Pairwise alignment for insertion detection

Each read is globally aligned to a target reference:
- Alignment blocks are parsed to detect where query advances more than target.
- These differences represent insertions in read relative to reference.
- We use score and percent identity, then positional checks near expected breakpoint for final validation.

### 2.5 MSA and weighted consensus

Per peak, insertion sequences are aligned together:
- Weighted input emphasizes recurrent sequences and reduces influence of single-read noise.
- Consensus at each column is called by weighted support.
- Columns with low support/coverage are called as `N`.

Why this matters:
- If two close but distinct ITDs are mixed into one peak, MSA columns become heterogeneous and consensus quality degrades (many `N`s). Local subpeak refinement is designed to prevent this.

### 2.6 Softmax normalization and z-score in validation

For each read, alignment-to-reference scores are transformed with softmax:
- Produces a per-read relative support profile across references.
- Higher `beta` sharpens winner-take-most behavior.

Then:
- Binary case (WT vs one ITD): direct probability/delta interpretation.
- Multi-reference case: z-score-like normalization of per-read support values before best-hit calling.

### 2.7 Strand bias testing

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

- **Merged close ITDs in one peak**: broad local mode leads to mixed MSA and ambiguous consensus.
- **Over-splitting peaks**: too-permissive splitting creates unstable low-read subpeaks.
- **Poor validation separation**: references too similar or low-quality reads reduce best-vs-second-best margin.
- **Strand imbalance artifacts**: primer/trimming or read-quality effects can mimic biology; compare with WT split.

## 5) Main Runtime Data Objects

- `reads_df`
  - Core per-read table: `read_id`, `read_seq`, `strand`, `read_len`, `gmm_peak_alias`, ambiguity flags.
- `comps`
  - Peak/component-level parameters: mean length, SD, fraction, aliases, AF-like metrics.
- `peak_subsets`
  - Mapping of peak alias to supporting read IDs.
- `all_itd_insertions` / insertion tables
  - Per-read insertion evidence: insertion position, length, sequence, alignment stats.
- `df_cons`
  - Per-peak consensus summary: consensus sequence and median insertion anchor.
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
