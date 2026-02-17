import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import argparse
import logging
import seaborn as sns
import shutil
import pymuscle5
import base64
import textwrap
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple
from venv import logger
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import Align, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sklearn.mixture import GaussianMixture
from scipy.stats import fisher_exact
from datetime import datetime


def export_itd_vcf(summary_df, ref_dict, itd_refs, comps, output_path):
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
    vcf_header = """##fileformat=VCFv4.2
##source=ITDValidator
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (ITD)">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##INFO=<ID=AF_GMM,Number=1,Type=Float,Description="Allele frequency from GMM clustering">
##INFO=<ID=AF_FITTED,Number=1,Type=Float,Description="Allele frequency from fitted model">
##INFO=<ID=AF_VALIDATED,Number=1,Type=Float,Description="Validated allele frequency">
##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Strand bias Fisher exact p-value">
##INFO=<ID=ITD_LEN,Number=1,Type=Integer,Description="Length of ITD insertion">
##INFO=<ID=INS_POS,Number=1,Type=Integer,Description="Insertion position on reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias p-value">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"""

    vcf_lines = [vcf_header]

    # --- Iterate over ITDs ---
    for _, row in summary_df.iterrows():
        itd = row["ref_alias"]
        if itd not in ref_dict or itd not in itd_refs:
            logger.warning(f"[export_itd_vcf] Skipping {itd}: missing reference entry.")
            continue

        ref_info = ref_dict[itd]
        ref_seq = str(ref_info["ref_seq_with_itd"])
        chr_name = ref_info.get("chr", "chrNA")
        pos = ref_info.get("genomic_insertion_pos", 0)
        alignment_pos = ref_info.get("local_insertion_pos", 0)
        ref_base = ref_seq[alignment_pos - 1] if 0 < alignment_pos <= len(ref_seq) else "N"
        alt_seq = itd_refs[itd].get("itd_seq", "N")
        itd_len = ref_info.get("itd_length", len(alt_seq))

        # --- Depth and allele frequencies ---
        dp = int(row.get("n_total_reads", 0))
        ad_alt = int(row.get("n_itd_reads", 0))
        ad_ref = max(dp - ad_alt, 0)
        af_val = float(row.get("allele_frequency", 0.0))
        sb_pval = row.get("fisher_p", np.nan)

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
            "AF_VALIDATED": round(af_val, 5),
            "SB_PVAL": round(sb_pval, 5) if not pd.isna(sb_pval) else ".",
            "ITD_LEN": itd_len,
            "INS_POS": pos
        }
        info_str = ";".join(f"{k}={v}" for k, v in info_fields.items())

        # --- FORMAT and SAMPLE ---
        gt = "0/1"
        format_str = "GT:DP:AF:AD:SB"
        sb_val = round(sb_pval, 5) if not pd.isna(sb_pval) else "."
        sample_str = f"{gt}:{dp}:{af_val:.4f}:{ad_ref},{ad_alt}:{sb_val}"

        vcf_line = f"{chr_name}\t{pos}\t{itd}\t{ref_base}\t{alt_seq}\t.\tPASS\t{info_str}\t{format_str}\t{sample_str}"
        vcf_lines.append(vcf_line)

    # --- Write to disk ---
    with open(output_path, "w") as f:
        f.write("\n".join(vcf_lines))

    logger.info(f"[export_itd_vcf] Wrote {len(summary_df)} ITDs to {output_path}")

def img_to_base64(path):
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        encoded = base64.b64encode(f.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded}"

def generate_itd_html_report(
    sample_name: str,
    reference_genome: str,
    seqio_reads: int,
    val_reads: int,
    summary_df: pd.DataFrame,
    itd_refs: dict,
    output_dir: str,
    plots_dir: str,
):
    """
    Generate an HTML summary report for ITD validation results.

    Parameters
    ----------
    sample_name : str
        Sample identifier
    reference_genome : str
        Reference genome build (e.g. "hg38")
    seqio_reads : int
        Total reads covering the region (from SeqIO)
    val_reads : int
        Number of validation-quality reads
    summary_df : pd.DataFrame
        DataFrame with ITD summary statistics (AF, DP, etc.)
    output_dir : str
        Output directory for HTML file
    plots_dir : str
        Directory where generated plots are stored
    """
    os.makedirs(output_dir, exist_ok=True)
    html_path = os.path.join(output_dir, f"{sample_name}_itd_report.html")

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    detected_itds = (len(summary_df["ref_alias"].tolist()) if len(summary_df) else "None")
    gmm_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_itd_gmm_fit_plot.png"))
    size_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_itd_size_distribution.png"))



    # HTML boilerplate
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>ITD Validation Report - {sample_name}</title>
<style>
body {{
  font-family: 'Arial', sans-serif;
  margin: 30px;
  background: #fafafa;
}}
h1, h2, h3 {{
  color: #2c3e50;
}}
.section {{
  background: #ffffff;
  padding: 15px 20px;
  margin-bottom: 25px;
  border-radius: 10px;
  box-shadow: 0px 2px 4px rgba(0,0,0,0.1);
}}
img {{
  max-width: 100%;
  border-radius: 6px;
  box-shadow: 0px 1px 3px rgba(0,0,0,0.2);
}}
.flex-row {{
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  gap: 20px;
}}
.plot-box {{
  flex: 1;
  min-width: 45%;
}}
</style>
</head>
<body>
<h1>ITD Validation Report</h1>
<div class="section">
  <h2>Sample Information</h2>
  <p><b>Sample name:</b> {sample_name}</p>
  <p><b>Reference genome:</b> {reference_genome}</p>
  <p><b>Generated on:</b> {timestamp}</p>
  <p><b>Total reads:</b> {seqio_reads:,}</p>
  <p><b>Good quality reads:</b> {val_reads:,}</p>
  <p><b>ITDs detected:</b> {detected_itds}</p>
</div>

<div class="section">
  <h2>General Overview</h2>
  <div class="flex-row">
    <div class="plot-box">
      <h3>Read Distribution and GMM Fit</h3>
      <img src="{gmm_plot}" alt="GMM fit plot">
    </div>
    <div class="plot-box">
      <h3>ITD Size Distribution</h3>
      <img src="{size_plot}" alt="ITD size distribution">
    </div>
  </div>
</div>
"""

    # --- Per ITD section ---
    for _, row in summary_df.iterrows():
        alias = row["ref_alias"]
        af = row.get("allele_frequency", 0)
        af = af * 100  # convert to percentage
        dp = row.get("n_total_reads", 0)
        itd_len = itd_refs.get(alias, {}).get("itd_length", "N/A")
        pos = itd_refs.get(alias, {}).get("genomic_insertion_pos", "N/A")
        msa_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_{alias}_MSA_consensus.png"))
        ref_plot = img_to_base64(os.path.join(plots_dir, f"{sample_name}_{alias}_itd_ref_hg38.png"))

        html += f"""
<div class="section">
  <h2>{alias}</h2>
  <table>
    <tr><th>Allele Frequency (Validated)</th><td>{af:.4f}%</td></tr>
    <tr><th>Depth (DP)</th><td>{dp}</td></tr>
    <tr><th>Insertion Length (bp)</th><td>{itd_len}</td></tr>
    <tr><th>Insertion Genomic Position</th><td>{pos}</td></tr>
  </table>
  <h3>Consensus Sequence Alignment</h3>
  <img src="{msa_plot}" alt="{alias} MSA plot">
  <h3>Genomic Context (hg38)</h3>
  <img src="{ref_plot}" alt="{alias} genomic reference">
</div>
"""

    html += """
</body>
</html>
"""

    with open(html_path, "w") as f:
        f.write(html)

    logger.info(f"[generate_itd_html_report] Wrote HTML report: {html_path}")
    return html_path

def call_no_itd(sample_name, genome, seqio_reads, output_folder, flt3_data_folder, temp_dir, html_report, logger, remove_intermediate_files):
            # Create empty VCF
        vcf_path = f"{sample_name}_FLT3_ITD_calls.vcf"
        output_path = os.path.join(output_folder, vcf_path)
        seqio_reads = len(seqio_reads)

        with open(output_path, "w") as vcf:
            vcf.write("""##fileformat=VCFv4.23
                        ##source=ITDValidator
                        ##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (ITD)">
                        ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">
                        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
                        ##INFO=<ID=AF_GMM,Number=1,Type=Float,Description="Allele frequency from GMM clustering">
                        ##INFO=<ID=AF_FITTED,Number=1,Type=Float,Description="Allele frequency from fitted model">
                        ##INFO=<ID=AF_VALIDATED,Number=1,Type=Float,Description="Validated allele frequency">
                        ##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Strand bias Fisher exact p-value">
                        ##INFO=<ID=ITD_LEN,Number=1,Type=Integer,Description="Length of ITD insertion">
                        ##INFO=<ID=INS_POS,Number=1,Type=Integer,Description="Insertion position on reference genome">
                        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
                        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">
                        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency (validated)">
                        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
                        ##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias p-value">
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE""")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            # no variants added
        logger.info(f"Empty VCF written: {output_path}")
        if html_report:
            html_name = f"{sample_name}_itd_report.html"
            output_path = os.path.join(output_folder, html_name)

            gmm_plot = img_to_base64(os.path.join(flt3_data_folder, f"{sample_name}_itd_gmm_fit_plot.png"))
            html = f"""
                        <!DOCTYPE html>
                        <html lang="en">
                        <head>
                        <meta charset="UTF-8">
                        <title>ITD Validation Report - {sample_name}</title>
                        <style>
                            body {{ font-family: Arial, sans-serif; margin: 40px; }}
                            .header {{ font-size: 24px; font-weight: bold; margin-bottom: 10px; }}
                            .subheader {{ font-size: 18px; margin-top: 20px; }}
                            .stats {{ font-family: monospace; margin-top: 10px; }}
                        </style>
                    </head>
                    <body>
                    <h1>ITD Validation Report</h1>
                    <div class="section">
                    <h2>Sample Information</h2>
                    <p><b>Sample name:</b> {sample_name}</p>
                    <p><b>Reference genome:</b> {genome}</p>
                    <p><b>Total reads:</b> {seqio_reads:,}</p>
                    <p><b>ITDs detected:</b> None </p>
                    </div>

                    <div class="section">
                    <h2>General Overview</h2>
                    <div class="flex-row">
                        <div class="plot-box">
                        <h3>Read Distribution and GMM Fit</h3>
                        <img src="{gmm_plot}" alt="GMM fit plot">
                        </div>
                    </div>
                    </div>
                """
            with open(output_path, "w") as f:
                f.write(html)
            logger.info(f"Empty HTML report generated: {output_path}")

        # Cleanup temp directory and intermediate files, if log is debug, keep all files
        if os.path.exists(temp_dir):
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug(f"Temp directory and FLT3 data folder retained for debugging: {temp_dir}")
            else:
                try:
                    shutil.rmtree(temp_dir)
                    logger.info(f"Cleaned up temp directory: {temp_dir}")
                except Exception as e:
                    logger.error(f"Error cleaning up temp directory: {e}")
        else:
            logger.warning(f"Temp directory does not exist, skipping cleanup: {temp_dir}")  
        
        if remove_intermediate_files:
            try:
                shutil.rmtree(flt3_data_folder)
                logger.info(f"Removed intermediate FLT3 data folder: {flt3_data_folder}")
            except Exception as e:
                logger.error(f"Error removing FLT3 data folder: {e}")