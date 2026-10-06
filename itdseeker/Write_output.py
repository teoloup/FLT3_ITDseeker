import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import numpy as np
import pandas as pd
import logging
import shutil
import base64
from datetime import datetime

logger = logging.getLogger(__name__)


# chr13 length per build, required for ##contig so bcftools can index the output.
CHR13_LENGTH = {"hg38": 114364328, "hg19": 115169878}


def _fmt_or(value):
    """Odds ratio for the VCF, guarding the degenerate all-one-strand case."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "."
    if v != v:            # NaN
        return "."
    if v in (float("inf"), float("-inf")):
        return "Inf"
    return f"{v:.3g}"


def strand_bias_is_real(odds_ratio, pval, plus, minus,
                        min_fold=3.0, max_p=0.05):
    """Decide whether a strand split is worth warning about.

    The Fisher p-value alone is not usable as a flag: it scales with depth, so a
    deep amplicon trips any fixed threshold on a difference far too small to
    matter. Measured on sample 10808, an ITD at 49.3% plus against a wild type at
    52.4% plus -- three percentage points, odds ratio 0.885 -- gives p = 0.0005 at
    13,537 reads and p = 0.71 at 269 reads. Same biology, opposite verdict.

    So the p-value is treated as necessary but not sufficient: the effect has to
    be a `min_fold` skew in the odds as well. A variant seen on essentially one
    strand is flagged regardless, since that is the pattern that actually
    indicates an artefact.
    """
    total = (plus or 0) + (minus or 0)
    if total:
        frac = (plus or 0) / total
        if frac <= 0.05 or frac >= 0.95:
            return True
    try:
        p = float(pval)
        orr = float(odds_ratio)
    except (TypeError, ValueError):
        return False
    if p != p or orr != orr or p > max_p:
        return False
    if orr in (float("inf"), float("-inf")) or orr == 0:
        return True
    return orr >= min_fold or orr <= 1.0 / min_fold


def format_pvalue(pval):
    """Render a p-value without collapsing small ones to 0.0."""
    if pval is None or pd.isna(pval):
        return "."
    return f"{float(pval):.3g}"


def build_vcf_header(genome_build=None, chr_name="chr13"):
    """Build the VCF header block shared by positive and negative calls.

    Assembled line by line rather than as a triple-quoted literal so the source
    file's own line endings can never leak stray CRs into the VCF.
    """
    lines = [
        "##fileformat=VCFv4.2",
        "##source=ITDValidator",
    ]
    if genome_build:
        lines.append(f"##reference={genome_build}")
        contig_len = CHR13_LENGTH.get(genome_build)
        if contig_len:
            lines.append(f"##contig=<ID={chr_name},length={contig_len}>")
    lines += [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (ITD)">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
        '##INFO=<ID=AF_GMM,Number=1,Type=Float,Description="Allele frequency from GMM clustering">',
        '##INFO=<ID=AF_FITTED,Number=1,Type=Float,Description="Allele frequency from fitted model">',
        '##INFO=<ID=FISHER_P,Number=1,Type=Float,Description="Fisher exact test p-value for strand bias vs the WT strand split. Depth-sensitive: interpret with STRAND_OR.">',
        '##INFO=<ID=STRAND_OR,Number=1,Type=Float,Description="Odds ratio of the ITD strand split against the WT split. 1.0 means no bias. This is the effect size that FISHER_P does not convey.">',
        '##INFO=<ID=ITD_LEN,Number=1,Type=Integer,Description="Length of ITD insertion">',
        '##INFO=<ID=INS_POS,Number=1,Type=Integer,Description="Insertion position on reference genome">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">',
        '##FORMAT=<ID=SB,Number=1,Type=String,Description="ITD-supporting strand counts as plus,minus">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    return "\n".join(lines)


def export_itd_vcf(summary_df, ref_dict, itd_refs, comps, output_path, genome_build=None):
    """
    Export detected ITDs into a standard VCF file with validated AF and strand bias info.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Table containing ITD summary (from calculate_allele_frequencies_and_strand_bias)
        Columns required: ref_alias, allele_frequency, n_itd_reads, n_total_reads, fisher_p
    ref_dict : dict
        {"ITD_1": {"chr": ..., "genomic_insertion_pos": ..., "ref_seq_with_itd": ..., "itd_length": ...}, ...}
    itd_refs : dict
        {"ITD_1": {"itd_seq": ..., "local_insertion_pos": ..., "genomic_insertion_pos": ..., "chr": ..., "itd_length": ...}, ...}
    comps : pd.DataFrame
        GMM component info with columns [peak_alias, effective_allele_freq, fraction]
    output_path : str
        Output VCF file path
    """

    logger.info(f"[export_itd_vcf] Exporting {len(summary_df)} ITDs to VCF: {output_path}")

    # --- VCF Header ---
    vcf_lines = [build_vcf_header(genome_build=genome_build)]

    # --- Iterate over ITDs ---
    for _, row in summary_df.iterrows():
        itd = row["ref_alias"]
        if itd not in ref_dict or itd not in itd_refs:
            logger.warning(f"[export_itd_vcf] Skipping {itd}: missing reference entry.")
            continue

        ref_info = ref_dict[itd]
        wt_seq = str(ref_dict.get("WT", {}).get("ref_seq_with_itd", ""))
        chr_name = ref_info.get("chr", "chrNA")
        pos_raw = itd_refs[itd].get("genomic_insertion_pos", ref_info.get("genomic_insertion_pos", 0))
        pos = int(pos_raw) if pd.notna(pos_raw) else 0
        alignment_pos = itd_refs[itd].get("local_insertion_pos", ref_info.get("local_insertion_pos"))
        if pd.notna(alignment_pos):
            alignment_pos = int(alignment_pos)
        ref_base = wt_seq[alignment_pos - 1] if wt_seq and isinstance(alignment_pos, int) and 0 < alignment_pos <= len(wt_seq) else "N"
        if ref_base == "N":
            logger.warning(
                f"[export_itd_vcf] REF base fallback to N for {itd}: "
                f"local_insertion_pos={alignment_pos}, wt_len={len(wt_seq)}, pos={pos}"
            )
        itd_seq = itd_refs[itd].get("itd_seq", "N")
        # VCF represents an insertion as REF=<anchor base>, ALT=<anchor base> +
        # <inserted bases>. Emitting the payload alone makes the record a 1bp->Nbp
        # substitution, which deletes the anchor base when a downstream tool
        # (bcftools norm, VEP) normalises it.
        alt_seq = f"{ref_base}{itd_seq}"
        itd_len = ref_info.get("itd_length", len(itd_seq))

        # --- Depth and allele frequencies ---
        dp = int(row.get("n_total_reads", 0))
        ad_alt = int(row.get("n_itd_reads", 0))
        ad_ref = max(dp - ad_alt, 0)
        af_val = float(row.get("allele_frequency", 0.0))
        sb_pval = row.get("fisher_p", np.nan)
        plus_raw = row.get("plus_reads", 0)
        minus_raw = row.get("minus_reads", 0)
        plus_reads = int(plus_raw) if pd.notna(plus_raw) else 0
        minus_reads = int(minus_raw) if pd.notna(minus_raw) else 0

        # --- GMM-derived frequencies ---
        af_gmm = comps.loc[comps["peak_alias"] == itd, "effective_allele_freq"].iloc[0] if itd in comps["peak_alias"].values else np.nan
        af_fit = comps.loc[comps["peak_alias"] == itd, "fraction"].iloc[0] if itd in comps["peak_alias"].values else np.nan

        # --- INFO ---
        info_fields = {
            "TYPE": "ITD",
            "AF": round(af_val, 5),
            "DP": dp,
            "AF_GMM": round(af_gmm, 5) if not pd.isna(af_gmm) else ".",
            "AF_FITTED": round(af_fit, 5) if not pd.isna(af_fit) else ".",
            "FISHER_P": format_pvalue(sb_pval),
            "STRAND_OR": _fmt_or(row.get("odds_ratio", None)),
            "ITD_LEN": itd_len,
            "INS_POS": pos
        }
        info_str = ";".join(f"{k}={v}" for k, v in info_fields.items())

        # --- FORMAT and SAMPLE ---
        gt = "0/1"
        format_str = "GT:DP:AF:AD:SB"
        sb_counts = f"{plus_reads},{minus_reads}"
        sample_str = f"{gt}:{dp}:{af_val:.4f}:{ad_ref},{ad_alt}:{sb_counts}"

        vcf_line = f"{chr_name}\t{pos}\t{itd}\t{ref_base}\t{alt_seq}\t.\tPASS\t{info_str}\t{format_str}\t{sample_str}"
        vcf_lines.append(vcf_line)

    # --- Write to disk ---
    # newline="\n" keeps the VCF byte-identical on Windows and Linux.
    with open(output_path, "w", newline="\n") as f:
        f.write("\n".join(vcf_lines) + "\n")

    logger.info(f"[export_itd_vcf] Wrote {len(summary_df)} ITDs to {output_path}")

def img_to_base64(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"

def _fmt_int(v):
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return str(v)


def _wrap_seq(seq, width=60):
    """Break a sequence into numbered lines, the way a sequence viewer would."""
    if not seq:
        return ""
    rows = []
    for i in range(0, len(seq), width):
        chunk = seq[i:i + width]
        # colour N separately: an N means the reads behind this consensus
        # disagreed at that column, which the reader needs to see
        marked = "".join(
            f'<span class="amb">{c}</span>' if c.upper() == "N" else c
            for c in chunk
        )
        rows.append(
            f'<div class="seqrow"><span class="seqpos">{i + 1:>5}</span>'
            f'<span class="seqbases">{marked}</span></div>'
        )
    return "".join(rows)


def _af_bar(af):
    pct = max(0.0, min(1.0, float(af))) * 100
    return (
        f'<div class="bar"><div class="barfill" style="width:{pct:.1f}%"></div></div>'
    )


def _strand_cell(plus, minus, pval, odds_ratio=None):
    total = (plus or 0) + (minus or 0)
    if not total:
        return '<span class="muted">n/a</span>'
    frac = (plus or 0) / total
    try:
        ptxt = f"{float(pval):.3g}"
    except (TypeError, ValueError):
        ptxt = "."
    ortxt = _fmt_or(odds_ratio)
    flag = ""
    if strand_bias_is_real(odds_ratio, pval, plus, minus):
        flag = ' <span class="warn">bias</span>'
    return (
        f'{_fmt_int(plus)} + / {_fmt_int(minus)} - '
        f'<span class="muted">({frac*100:.1f}% plus)</span>'
        f'<div class="ministack"><div class="ministack-plus" style="width:{frac*100:.0f}%"></div></div>'
        f'<span class="muted">odds ratio {ortxt}, Fisher p={ptxt}</span>{flag}'
    )


def generate_itd_html_report(
    sample_name: str,
    reference_genome: str,
    seqio_reads: int,
    val_reads: int,
    summary_df: pd.DataFrame,
    itd_refs: dict,
    output_dir: str,
    plots_dir: str,
    df_cons: pd.DataFrame = None,
    haplotype_method: str = None,
    rescue_note: str = None,
):
    """Write a self-contained HTML report for one sample.

    Everything is inlined (images as data URIs) so the file can be emailed or
    archived on its own. Reports one card per ITD carrying the numbers a reader
    actually needs to act on: size, position, allele frequency, read support,
    strand balance, and the inserted sequence itself.
    """
    os.makedirs(output_dir, exist_ok=True)
    html_path = os.path.join(output_dir, f"{sample_name}_itd_report.html")
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    n_itds = len(summary_df) if summary_df is not None else 0
    gmm_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_itd_gmm_fit_plot.png"))
    size_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_itd_size_distribution.png"))

    cons_by_alias = {}
    if df_cons is not None and not df_cons.empty:
        for _, r in df_cons.iterrows():
            cons_by_alias[str(r["peak_alias"])] = r

    total_af = 0.0
    if n_itds:
        try:
            total_af = float(summary_df["allele_frequency"].sum())
        except Exception:
            total_af = 0.0

    verdict = "POSITIVE" if n_itds else "NEGATIVE"
    verdict_class = "pos" if n_itds else "neg"

    # ---------------------------------------------------------------- cards --
    cards = []
    for _, row in (summary_df.iterrows() if n_itds else []):
        alias = str(row["ref_alias"])
        info = itd_refs.get(alias, {})
        itd_seq = info.get("itd_seq", "") or ""
        itd_len = info.get("itd_length", len(itd_seq))
        gpos = info.get("genomic_insertion_pos", "n/a")
        chrom = info.get("chr", "chr13")
        af = float(row.get("allele_frequency", 0.0))
        dp = int(row.get("n_total_reads", 0) or 0)
        alt_n = int(row.get("n_itd_reads", 0) or 0)
        ref_n = max(dp - alt_n, 0)
        plus = row.get("plus_reads", 0)
        minus = row.get("minus_reads", 0)
        pval = row.get("fisher_p", None)
        orr = row.get("odds_ratio", None)

        cons = cons_by_alias.get(alias)
        n_count = int(cons["consensus_seq"].count("N")) if cons is not None else itd_seq.count("N")
        amb_note = ""
        if n_count:
            amb_note = (
                f'<div class="note warnbox"><b>{n_count} ambiguous base'
                f'{"" if n_count == 1 else "s"} (N)</b> in this consensus. The reads '
                f'behind this peak disagree, which usually means more than one ITD of '
                f'the same length is present and has not been separated.</div>'
            )

        flags = []
        if n_count:
            flags.append(('warn', f'{n_count} ambiguous base' + ('' if n_count == 1 else 's')))
        if dp and alt_n < 50:
            flags.append(('warn', f'only {alt_n} supporting reads'))
        if strand_bias_is_real(orr, pval, plus, minus):
            flags.append(('warn', 'strand bias'))
        if af < 0.02:
            flags.append(('info', 'low allele frequency'))
        if itd_len % 3:
            flags.append(('info', 'frameshift'))
        flag_html = "".join(
            f'<span class="flag {k}">{t}</span>' for k, t in flags
        ) or '<span class="flag ok">no flags</span>'

        in_frame = (itd_len % 3 == 0)
        frame_txt = ("in frame" if in_frame else
                     f'<span class="warn">out of frame ({itd_len % 3})</span>')

        msa_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_{alias}_MSA_consensus.png"))
        ref_plot = img_to_base64(
            os.path.join(plots_dir, f"{sample_name}_{alias}_itd_ref_{reference_genome}.png")
        )
        pileup = img_to_base64(os.path.join(plots_dir, f"{sample_name}_{alias}_pileup.png"))
        plots = ""
        if pileup:
            plots += (f'<figure><figcaption>Supporting reads, one row each, with the '
                      f'inserted bases marked and coloured by strand. A clean vertical '
                      f'edge means every read placed the insertion at the same position.'
                      f'</figcaption><img src="{pileup}" alt="{alias} read pileup"></figure>')
        if ref_plot:
            plots += (f'<figure><figcaption>Position within FLT3</figcaption>'
                      f'<img src="{ref_plot}" alt="{alias} genomic context"></figure>')
        if msa_plot:
            plots += (f'<figure><figcaption>Consensus support per alignment column</figcaption>'
                      f'<img src="{msa_plot}" alt="{alias} consensus coverage"></figure>')

        cards.append(f"""
<section class="card">
  <header class="cardhead">
    <h3>{alias}</h3>
    <div class="chips">
      <span class="chip">{itd_len} bp</span>
      <span class="chip">{frame_txt}</span>
      <span class="chip mono">{chrom}:{_fmt_int(gpos)}</span>
    </div>
  </header>
  <div class="flags">{flag_html}</div>

  <div class="grid">
    <div class="metric">
      <div class="label">Allele frequency</div>
      <div class="value big">{af * 100:.2f}<span class="unit">%</span></div>
      {_af_bar(af)}
      <div class="muted">{_fmt_int(alt_n)} of {_fmt_int(dp)} classified reads</div>
    </div>
    <div class="metric">
      <div class="label">Read support</div>
      <div class="value">{_fmt_int(alt_n)}<span class="unit"> ITD</span></div>
      <div class="muted">{_fmt_int(ref_n)} reference-supporting</div>
      <div class="muted">depth {_fmt_int(dp)}</div>
    </div>
    <div class="metric">
      <div class="label">Strand balance</div>
      <div class="value small">{_strand_cell(plus, minus, pval, orr)}</div>
    </div>
  </div>

  {amb_note}

  <details open>
    <summary>Inserted sequence &mdash; {itd_len} bp</summary>
    <div class="note">
      Shown on the <b>plus strand of {chrom}</b>, matching the VCF. FLT3 is
      transcribed from the minus strand, so this is the reverse complement of the
      coding sequence.
    </div>
    <div class="seqbox">{_wrap_seq(itd_seq)}</div>
    <div class="muted">
      VCF ALT is this sequence prefixed by the reference base at
      {chrom}:{_fmt_int(gpos)}.
    </div>
  </details>

  {f'<div class="figures">{plots}</div>' if plots else ''}
</section>""")

    if not n_itds:
        cards.append("""
<section class="card">
  <h3>No ITD reported</h3>
  <p class="muted">
    No insertion passed validation and the allele-frequency threshold. Check the
    read-length distribution below: a sample with too few reads, or no peak above
    the wild-type amplicon length, cannot yield a call.
  </p>
</section>""")

    method_line = ""
    if haplotype_method:
        method_line = f'<div><dt>Haplotype method</dt><dd class="mono">{haplotype_method}</dd></div>'
    rescue_line = ""
    if rescue_note:
        rescue_line = f'<div class="note warnbox">{rescue_note}</div>'

    overview = ""
    if gmm_plot:
        overview += (f'<figure><figcaption>Read-length distribution and fitted peaks</figcaption>'
                     f'<img src="{gmm_plot}" alt="GMM fit"></figure>')
    if size_plot:
        overview += (f'<figure><figcaption>Insertion sizes across reads</figcaption>'
                     f'<img src="{size_plot}" alt="ITD size distribution"></figure>')

    used_pct = (100.0 * val_reads / seqio_reads) if seqio_reads else 0.0

    # One row per ITD, so a reader can take in the whole sample before reading
    # any individual card.
    summary_table = ""
    if n_itds:
        trs = []
        for _, row in summary_df.iterrows():
            a = str(row["ref_alias"])
            inf = itd_refs.get(a, {})
            ln = inf.get("itd_length", 0)
            trs.append(
                f'<tr><td class="mono">{a}</td>'
                f'<td>{ln} bp</td>'
                f'<td>{"in frame" if ln % 3 == 0 else "frameshift"}</td>'
                f'<td class="mono">{inf.get("chr","chr13")}:{_fmt_int(inf.get("genomic_insertion_pos","n/a"))}</td>'
                f'<td class="num"><b>{float(row.get("allele_frequency",0))*100:.2f}%</b></td>'
                f'<td class="num">{_fmt_int(row.get("n_itd_reads",0))}</td>'
                f'<td class="num">{_fmt_int(row.get("n_total_reads",0))}</td></tr>'
            )
        summary_table = f"""
  <h2>Summary</h2>
  <section class="card">
    <table class="summary">
      <thead><tr><th>ITD</th><th>Size</th><th>Frame</th><th>Position</th>
        <th class="num">AF</th><th class="num">Reads</th><th class="num">Depth</th></tr></thead>
      <tbody>{''.join(trs)}</tbody>
    </table>
  </section>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>FLT3-ITD report &mdash; {sample_name}</title>
<style>
  :root {{
    --bg:#f6f7f9; --panel:#fff; --ink:#1c2430; --muted:#6b7684;
    --line:#e3e7ec; --accent:#2f6f4f; --accent-soft:#e7f2ec;
    --warn:#9a4b17; --warn-soft:#fdf1e6; --neg:#54606e;
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; padding:24px 16px 56px; background:var(--bg); color:var(--ink);
    font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
  .wrap {{ max-width:1040px; margin:0 auto; }}
  h1 {{ font-size:22px; margin:0 0 2px; }}
  h2 {{ font-size:15px; text-transform:uppercase; letter-spacing:.07em;
    color:var(--muted); margin:32px 0 12px; font-weight:600; }}
  h3 {{ font-size:17px; margin:0; }}
  .mono {{ font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; }}
  .muted {{ color:var(--muted); font-size:13px; }}
  .warn {{ color:var(--warn); font-weight:600; }}

  .masthead {{ display:flex; justify-content:space-between; align-items:flex-start;
    gap:16px; flex-wrap:wrap; border-bottom:2px solid var(--ink); padding-bottom:14px; }}
  .verdict {{ font-size:13px; font-weight:700; letter-spacing:.08em;
    padding:6px 14px; border-radius:999px; white-space:nowrap; }}
  .verdict.pos {{ background:var(--accent-soft); color:var(--accent);
    border:1px solid var(--accent); }}
  .verdict.neg {{ background:#eef1f4; color:var(--neg); border:1px solid var(--neg); }}

  dl.meta {{ display:flex; flex-wrap:wrap; gap:0 32px; margin:16px 0 0; }}
  dl.meta div {{ min-width:120px; }}
  dt {{ font-size:12px; text-transform:uppercase; letter-spacing:.05em;
    color:var(--muted); }}
  dd {{ margin:2px 0 10px; font-weight:600; }}

  .card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
    padding:20px 22px; margin-bottom:18px; }}
  .cardhead {{ display:flex; align-items:center; gap:14px; flex-wrap:wrap;
    border-bottom:1px solid var(--line); padding-bottom:12px; margin-bottom:16px; }}
  .chips {{ display:flex; gap:8px; flex-wrap:wrap; }}
  .chip {{ font-size:12px; background:#eef1f4; border-radius:5px; padding:3px 9px; }}
  .chip.mono {{ font-family:ui-monospace,Menlo,Consolas,monospace; }}

  .grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr));
    gap:20px; }}
  .label {{ font-size:12px; text-transform:uppercase; letter-spacing:.05em;
    color:var(--muted); margin-bottom:4px; }}
  .value {{ font-size:20px; font-weight:650; }}
  .value.big {{ font-size:30px; line-height:1.1; }}
  .value.small {{ font-size:14px; font-weight:500; }}
  .unit {{ font-size:14px; font-weight:500; color:var(--muted); }}

  .bar {{ height:7px; background:#e7eaee; border-radius:4px; overflow:hidden;
    margin:8px 0 6px; }}
  .barfill {{ height:100%; background:var(--accent); }}
  .ministack {{ height:5px; background:#d8626f; border-radius:3px; overflow:hidden;
    margin:6px 0 4px; max-width:160px; }}
  .ministack-plus {{ height:100%; background:#4a7fb5; }}

  details {{ margin-top:18px; border-top:1px solid var(--line); padding-top:14px; }}
  summary {{ cursor:pointer; font-weight:600; font-size:14px; }}
  .note {{ font-size:13px; color:var(--muted); margin:10px 0; }}
  .warnbox {{ background:var(--warn-soft); border-left:3px solid var(--warn);
    color:var(--warn); padding:10px 12px; border-radius:0 6px 6px 0; }}

  .seqbox {{ background:#fbfcfd; border:1px solid var(--line); border-radius:6px;
    padding:12px 14px; overflow-x:auto; margin:10px 0; }}
  .seqrow {{ font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
    font-size:13px; white-space:pre; letter-spacing:.06em; }}
  .seqpos {{ color:var(--muted); margin-right:14px; user-select:none; }}
  .seqbases {{ word-break:break-all; }}
  .amb {{ background:var(--warn-soft); color:var(--warn); font-weight:700; }}

  .flags {{ display:flex; gap:6px; flex-wrap:wrap; margin:-6px 0 14px; }}
  .flag {{ font-size:11px; padding:2px 8px; border-radius:999px; font-weight:600; }}
  .flag.warn {{ background:var(--warn-soft); color:var(--warn); }}
  .flag.info {{ background:#eef1f4; color:var(--muted); }}
  .flag.ok {{ background:var(--accent-soft); color:var(--accent); }}

  table.summary {{ width:100%; border-collapse:collapse; font-size:14px; }}
  table.summary th {{ text-align:left; font-size:11px; text-transform:uppercase;
    letter-spacing:.05em; color:var(--muted); border-bottom:1px solid var(--line);
    padding:0 10px 7px 0; font-weight:600; }}
  table.summary td {{ padding:9px 10px 9px 0; border-bottom:1px solid var(--line); }}
  table.summary tr:last-child td {{ border-bottom:none; }}
  .num {{ text-align:right; }}

  .figures {{ display:grid; grid-template-columns:1fr; gap:18px; margin-top:18px; }}
  figure {{ margin:0; }}
  figcaption {{ font-size:12px; color:var(--muted); margin-bottom:6px; }}
  img {{ max-width:100%; border:1px solid var(--line); border-radius:6px; display:block; }}

  @media print {{
    body {{ background:#fff; padding:0; }}
    .card {{ break-inside:avoid; border-color:#ccc; }}
    details {{ display:block; }}
  }}
</style>
</head>
<body>
<div class="wrap">

  <div class="masthead">
    <div>
      <h1>FLT3-ITD report</h1>
      <div class="muted">{sample_name} &middot; {reference_genome} &middot; {timestamp}</div>
    </div>
    <div class="verdict {verdict_class}">{verdict}
      {f'&middot; {n_itds} ITD' + ('' if n_itds == 1 else 's') if n_itds else ''}</div>
  </div>

  <dl class="meta">
    <div><dt>Reads in region</dt><dd>{_fmt_int(seqio_reads)}</dd></div>
    <div><dt>Passed to validation</dt><dd>{_fmt_int(val_reads)} <span class="muted">({used_pct:.1f}%)</span></dd></div>
    <div><dt>ITDs reported</dt><dd>{n_itds}</dd></div>
    <div><dt>Combined ITD burden</dt><dd>{total_af * 100:.2f}%</dd></div>
    {method_line}
  </dl>

  {rescue_line}

  {summary_table}

  <h2>Detected ITDs</h2>
  {''.join(cards)}

  <h2>Supporting evidence</h2>
  <section class="card">
    <div class="figures">{overview or '<p class="muted">No overview plots were produced.</p>'}</div>
  </section>

  <p class="muted" style="margin-top:28px">
    Allele frequency is the share of validated, classified reads assigned to this
    ITD by competitive alignment against the wild-type and per-ITD references.
    Strand balance compares this ITD plus/minus split against the wild-type split.
    The odds ratio is the effect size (1.0 means no bias); the Fisher p-value
    scales with depth, so on a deep amplicon a three-percentage-point difference
    can reach p&lt;0.001 while being biologically meaningless. A bias flag needs
    both: a significant p AND at least a three-fold skew in the odds, or a
    variant seen on essentially one strand. Positions are 1-based on the plus strand.
  </p>
</div>
</body>
</html>
"""

    with open(html_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(html)

    logger.info(f"[generate_itd_html_report] Wrote HTML report: {html_path}")
    return html_path


def call_no_itd(
    sample_name,
    genome,
    seqio_reads,
    output_folder,
    flt3_data_folder,
    temp_dir,
    html_report,
    logger=None,
    remove_intermediate_files=False,
    reason="No ITDs detected.",
    temp_dir_is_ours=True,
):
    if logger is None:
        logger = logging.getLogger(__name__)

    vcf_path = f"{sample_name}_FLT3_ITD_calls.vcf"
    output_path = os.path.join(output_folder, vcf_path)
    total_reads = seqio_reads if isinstance(seqio_reads, int) else len(seqio_reads)

    with open(output_path, "w", newline="\n") as vcf:
        vcf.write(build_vcf_header(genome_build=genome) + "\n")
    logger.info(f"Empty VCF written: {output_path}")

    if html_report:
        html_name = f"{sample_name}_itd_report.html"
        output_path = os.path.join(output_folder, html_name)

        gmm_plot = img_to_base64(os.path.join(flt3_data_folder, f"{sample_name}_itd_gmm_fit_plot.png"))
        gmm_plot_block = (
            f'<img src="{gmm_plot}" alt="GMM fit plot">'
            if gmm_plot
            else "<p><i>Read-length / GMM plot not available for this negative run.</i></p>"
        )
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>ITD Validation Report - {sample_name}</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 40px; background: #fafafa; }}
.section {{ background: #ffffff; padding: 15px 20px; margin-bottom: 25px; border-radius: 10px; box-shadow: 0px 2px 4px rgba(0,0,0,0.1); }}
img {{ max-width: 100%; border-radius: 6px; box-shadow: 0px 1px 3px rgba(0,0,0,0.2); }}
</style>
</head>
<body>
<h1>ITD Validation Report</h1>
<div class="section">
  <h2>Sample Information</h2>
  <p><b>Sample name:</b> {sample_name}</p>
  <p><b>Reference genome:</b> {genome}</p>
  <p><b>Total reads:</b> {total_reads:,}</p>
  <p><b>ITDs detected:</b> 0</p>
  <p><b>Result:</b> No ITDs detected.</p>
  <p><b>Reason:</b> {reason}</p>
</div>
<div class="section">
  <h2>General Overview</h2>
  <h3>Read Distribution and GMM Fit</h3>
  {gmm_plot_block}
</div>
</body>
</html>
"""
        with open(output_path, "w") as f:
            f.write(html)
        logger.info(f"Empty HTML report generated: {output_path}")

    if os.path.exists(temp_dir) and temp_dir_is_ours:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Temp directory and FLT3 data folder retained for debugging: {temp_dir}")
        else:
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temp directory: {temp_dir}")
            except Exception as e:
                logger.error(f"Error cleaning up temp directory: {e}")
    else:
        logger.debug(f"Leaving temp directory in place: {temp_dir}")

    if remove_intermediate_files:
        try:
            shutil.rmtree(flt3_data_folder)
            logger.info(f"Removed intermediate FLT3 data folder: {flt3_data_folder}")
        except Exception as e:
            logger.error(f"Error removing FLT3 data folder: {e}")
