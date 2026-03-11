import os
os.environ["MPLBACKEND"] = "Agg"       # disable any GUI backend
os.environ["DISPLAY"] = ""             # make sure Tk can't open a window
os.environ["TK_SILENCE_DEPRECATION"] = "1"
import pandas as pd
import argparse
import logging
import shutil
import sys
import traceback
from pathlib import Path
from Bio.Seq import Seq

from bam_extractor import extract_flt3_reads
from GMM_peaks import fit_gmm_itds, plot_gmm_itds, refine_peak_substructure_once
from Pairwise_aligment_toolkit import align_reads_multi_ref_parallel
from Write_output import export_itd_vcf, generate_itd_html_report, call_no_itd
from Helper_functions import extract_itd_insertions_from_subset_parallel, plot_itd_size_distribution, build_itd_reference_per_peak, make_validation_refs, prepare_validation_reads, calculate_allele_frequencies_and_strand_bias
from Multiple_seq_aligment_toolkit import build_itd_consensus_sequences

def install_unhandled_exception_logger():
    """
    Install a global exception hook so unexpected failures are logged clearly.
    This is for true error conditions and does not emit negative/no-ITD outputs.
    """
    def _hook(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        logger = logging.getLogger()
        tb = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback)).rstrip()
        logger.critical("=" * 90)
        logger.critical("FATAL: Unhandled exception. Aborting run.")
        logger.critical("No negative VCF/report is written for this error condition.")
        logger.critical(f"Exception: {exc_type.__name__}: {exc_value}")
        logger.critical("Traceback follows:\n%s", tb)
        logger.critical("=" * 90)

    sys.excepthook = _hook

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect and validate FLT3-ITD events from Nanopore amplicon BAM files."
    )
    parser.add_argument(
        "-b", "--bam", type=str, required=True, help="Path to input BAM file."   
    )
    parser.add_argument(
        "-o", "--output-folder", type=str, required=True, help="Output folder."
    )
    parser.add_argument(
        "--html-report", action="store_true", help="Create also an html report."
    )
    parser.add_argument(
        "-s", "--sample-name", type=str, required=True, help="Sample name to add as prefix to output files."
    )
    parser.add_argument(
        "-g", "--genome", type=str, default="hg38", choices=['hg38','hg19'], help="Human genome version the bam is aligned"
    )
    parser.add_argument(
        "--min-allele-frequency", type=float, default=0.05, help="Minimum allele frequency of ITD to report (default: 0.05)."
    )
    parser.add_argument(
        "-q", "--min-mapping-quality",  type = int, default = 20, help = "Minimum mapping quality of reads to be extracted from the FLT3 region of the BAM file (default: 20)"
    )
    parser.add_argument(
        "--min-read-length",  type = int, default = 330, help = "Minimum read length of reads after primer trimming (current amplicon length of wild type 336)"
    )
    parser.add_argument(
        "--wt-amplicon-length",  type = int, default = 336, help = "Amplicon length of wild type"
    )
    parser.add_argument(
        "--per-peak-read-assignment-mode",  type = str, default = "manual", choices=['manual', 'predict_proba', 'hybrid'], help = "Read assignment mode for each peak (default: manual)"
    )
    parser.add_argument(
        "--probab-threshold-for-peak-assignment",  type = float, default = 0.85, help = "Probability threshold for peak assignment when using predict_proba or hybrid mode (default: 0.85)"
    )
    parser.add_argument(
        "--max-itd-size", type=int, default=300, help="Maximum ITD size (default: 300).",
    )
    parser.add_argument(
        "--min-itd-size", type=int, default=12, help="Minimum ITD size (default: 12).",
    )
    parser.add_argument(
        "--max-itds-to-detect", type=int, default=6, help="Maximum number of ITD to validate and report (default: 6).",
    )
    parser.add_argument(
        "--min-gmm-fraction", type=float, default=0.01, help="Minimum allele frequency to keep a GMM component (default: 0.01).",
    )
    parser.add_argument(
        "--min-gmm-peak-distance", type=int, default=10, help="Minimum distance between GMM peaks (default: 10).",
    )
    parser.add_argument(
        "--disable-subpeak-refinement", action="store_true", help="Disable one-level local refinement of each initial ITD peak."
    )
    parser.add_argument(
        "--min-reads-for-subpeak-refinement", type=int, default=150, help="Minimum reads in a parent peak to attempt one-level subpeak refinement (default: 150)."
    )
    parser.add_argument(
        "--min-subpeak-fraction", type=float, default=0.15, help="Minimum fraction per child subpeak when splitting a parent peak (default: 0.15)."
    )
    parser.add_argument(
        "--min-subpeak-distance", type=float, default=3.0, help="Minimum distance (bp) between child means to keep a split (default: 3.0)."
    )
    parser.add_argument(
        "--max-subpeak-sd", type=float, default=5.0, help="Maximum SD (bp) for each child subpeak in refinement (default: 5.0)."
    )
    parser.add_argument(
        "--min-bic-gain-for-subpeak-split", type=float, default=10.0, help="Minimum BIC gain (k=1 minus k=2) to accept a local split (default: 10.0)."
    )
    parser.add_argument(
        "--msa-max-unique", type=int, default=300, help="Maximum unique insertion sequences used per peak for MSA consensus (default: 300)."
    )
    parser.add_argument(
        "--msa-min-weight-coverage", type=float, default=0.98, help="Minimum cumulative weight coverage when selecting insertion panel for MSA (default: 0.98)."
    )
    parser.add_argument(
        "--msa-base-threshold", type=float, default=0.7, help="Minimum weighted base fraction to emit a non-ambiguous base in consensus (default: 0.7)."
    )
    parser.add_argument(
        "--msa-min-col-coverage", type=float, default=0.7, help="Minimum weighted non-gap coverage required to keep an MSA column in consensus (default: 0.7)."
    )
    parser.add_argument(
        "--max-peak-sd", type=float, default=5.0, help="Maximum standard deviation for GMM peaks (default: 5.0).",
    )
    parser.add_argument(
        "--force-number-of-peaks", type=int, default=None, help="If set, force this number of components instead of using BIC.",
    )
    parser.add_argument(
        "--temp-dir", type=str, default=None, help="Path to create a temp dir for intermediate files, (default: the output-folder)",
    )
    parser.add_argument(
        "--remove-intermediate-files", action="store_true", help="If set, remove intermediate files after run."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1, help="Number of threads to use (default: 1)."
    )
    parser.add_argument(
        "--log-level", type=str, default="INFO", choices=['INFO','DEBUG'], help="Logging level (default: INFO)."
    )

    args = parser.parse_args()
    
    # creating the logger object, and setting the log level
    logger = logging.getLogger()
    logger.setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    if not logger.hasHandlers():  # prevents double logging if run multiple times
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        ))
        logger.addHandler(handler)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    install_unhandled_exception_logger()
    
    
    # Read args and declare variables
    bam_file = args.bam
    output_folder = Path(args.output_folder)
    flt3_data_folder = os.path.join(output_folder, "flt3_data")
    flt3_data_folder = Path(flt3_data_folder)
    html_report = args.html_report
    sample_name = args.sample_name
    min_itd_size = args.min_itd_size
    max_itds_detected = args.max_itds_to_detect
    min_allele_frequency = args.min_allele_frequency
    min_gmm_fraction = args.min_gmm_fraction
    max_peak_sd = args.max_peak_sd
    min_gmm_peak_distance = args.min_gmm_peak_distance
    enable_subpeak_refinement = not args.disable_subpeak_refinement
    min_reads_for_subpeak_refinement = args.min_reads_for_subpeak_refinement
    min_subpeak_fraction = args.min_subpeak_fraction
    min_subpeak_distance = args.min_subpeak_distance
    max_subpeak_sd = args.max_subpeak_sd
    min_bic_gain_for_subpeak_split = args.min_bic_gain_for_subpeak_split
    msa_max_unique = args.msa_max_unique
    msa_min_weight_coverage = args.msa_min_weight_coverage
    msa_base_threshold = args.msa_base_threshold
    msa_min_col_coverage = args.msa_min_col_coverage
    min_mapq = args.min_mapping_quality
    min_read_length = args.min_read_length
    wt_amplicon_length = args.wt_amplicon_length
    genome=args.genome
    force_k = args.force_number_of_peaks
    threads = args.threads
    max_itd_length = args.max_itd_size
    temp_dir = args.temp_dir
    peak_read_assignment_mode = args.per_peak_read_assignment_mode
    prob_threshold = args.probab_threshold_for_peak_assignment
    remove_intermediate_files = args.remove_intermediate_files

    #basic input checks, folder creations etc
    if not output_folder.exists():
        try:
            logger.info(f"Creating output folder: {output_folder}")
            output_folder.mkdir(parents=True, exist_ok=True)
            logger.info(f"Output folder: {output_folder}")
        except Exception as e:
            logger.error(f"Error creating output folder: {e}")
            quit(1)
    if not flt3_data_folder.exists():
        try:
            logger.info(f"Creating FLT3 data folder: {flt3_data_folder}")
            flt3_data_folder.mkdir(parents=True, exist_ok=True)
            logger.info(f"FLT3 data folder: {flt3_data_folder}")
        except Exception as e:
            logger.error(f"Error creating FLT3 data folder: {e}")
            quit(1)

    if not Path(bam_file).exists():
        logger.error(f"BAM file does not exist: {bam_file}")
        quit(1)

    if temp_dir is None:
        temp_dir = os.path.join(output_folder, "temp")

    logger.info(f"Sample name: {sample_name}")
    logger.info(f"Human genome assembly: {genome}")
    logger.info(f"Min ITD size: {min_itd_size}")
    logger.info(f"Max ITD size: {max_itd_length}")
    logger.info(f"Max ITDs to detect: {max_itds_detected}")
    logger.info(f"Min GMM fraction: {min_gmm_fraction}")
    logger.info(f"Min allele frequency of ITD: {min_allele_frequency}")
    logger.info(f"Force number of peaks: {force_k}")
    logger.info(f"Max STD of peaks: {max_peak_sd}")
    logger.info(f"Method to assign reads to peaks: {peak_read_assignment_mode}")
    logger.info(f"Subpeak refinement enabled: {enable_subpeak_refinement}")
    logger.info(f"MSA max unique panel size: {msa_max_unique}")
    logger.info(f"Number of threads requested: {threads}")
    logger.debug(f"Temp directory path: {temp_dir}")

    def finalize_no_itd(reason):
        logger.warning(reason)
        call_no_itd(
            sample_name,
            genome,
            seqio_reads,
            output_folder,
            flt3_data_folder,
            temp_dir,
            html_report,
            logger,
            remove_intermediate_files,
            reason=reason,
        )
        logger.info("Exiting program with empty outputs (no ITDs detected).")
        raise SystemExit(0)

    # FLT3 reference amplicon and genomic coordinates used by the pipeline.
    DEFAULT_REF_WT = """CTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCT"""
    ref_seq = Seq(DEFAULT_REF_WT)

    amplicon_coords = {
    "hg19": {"chr": "chr13", "start": 28608018, "end": 28608353},
    "hg38": {"chr": "chr13", "start": 28033881, "end": 28034216},
    }
    flt3_exons_hg38 = [
    (28003273, 28004174),(28014451, 28014557),(28015156, 28015256),(28015589, 28015701),(28018466, 28018589),(28023349, 28023477),(28024860, 28024943),(28027087, 28027241),
    (28028177, 28028288),(28033886, 28033991),(28034081, 28034214),(28034300, 28034407),(28035494, 28035673),(28035934, 28036043),(28037184, 28037288),(28048274, 28048443),
    (28049383, 28049537),(28049634, 28049774),(28050094, 28050222),(28052544, 28052674),(28057346, 28057462),(28061866, 28062069),(28070490, 28070612),(28100467, 28100576)
    ]

    flt3_exons_hg19 = [
    (28577410, 28578311),(28588588, 28588694),(28589293, 28589393),(28589726, 28589838),(28592603, 28592726),(28597486, 28597614),(28598997, 28599080),(28601224, 28601378),
    (28602314, 28602425),(28608023, 28608128),(28608218, 28608351),(28608437, 28608544),(28609631, 28609810),(28610071, 28610180),(28611321, 28611425),(28622411, 28622580),
    (28623520, 28623674),(28623771, 28623911),(28624231, 28624359),(28626681, 28626811),(28631483, 28631599),(28636003, 28636206),(28644627, 28644749),(28674604, 28674713)
    ]

    if genome == 'hg38':
        exon_boundaries = flt3_exons_hg38
    elif genome == 'hg19':
        exon_boundaries = flt3_exons_hg19
    else:
        logger.error(f"Unsupported genome build: {genome}")
        quit(1)
    exon_labels = [f"Ex{idx+1}" for idx in range(len(exon_boundaries))]

    # Extract FLT3 reads
    try:
        seqio_reads , trimmed_fastq = extract_flt3_reads(
            bam_file=bam_file,
            genome_build=genome,
            threads=threads,
            min_mapping_quality = min_mapq,
            min_length = min_read_length,
            max_itd_length = max_itd_length,
            wt_amplicon_length = wt_amplicon_length,
            temp_dir = temp_dir
        )
    except FileNotFoundError as e:
        logger.error(f"Required external tool not found: {e}")
        logger.error("Aborting without writing negative (no-ITD) outputs.")
        quit(1)
    except Exception as e:
        logger.error(f"Read extraction/trimming failed: {e}")
        logger.error("Aborting without writing negative (no-ITD) outputs.")
        quit(1)

    if not seqio_reads:
        finalize_no_itd("No FLT3 reads available after extraction/trimming.")

    try:
        gmm_fit = fit_gmm_itds(
            seqio_reads,
            assign_mode=peak_read_assignment_mode,
            prob_threshold=prob_threshold,
            assign_width_factor=2.0,
            min_gmm_fraction=min_gmm_fraction,
            max_peak_sd=max_peak_sd,
            max_itds_detected=max_itds_detected,
            min_ggmm_peak_distance=min_gmm_peak_distance,
            force_k=force_k
        )
    except RuntimeError as e:
        if "No GMM components passed filtering criteria" in str(e):
            finalize_no_itd(str(e))
        logger.error(f"GMM fitting failed: {e}")
        logger.error("Aborting without writing negative (no-ITD) outputs.")
        quit(1)
    except Exception as e:
        logger.error(f"GMM fitting failed: {e}")
        logger.error("Aborting without writing negative (no-ITD) outputs.")
        quit(1)
    gmm = gmm_fit.gmm
    comps = gmm_fit.comps
    reads_df = gmm_fit.reads_df
    peak_subsets = gmm_fit.peak_subsets

    if enable_subpeak_refinement:
        logger.info("Running one-level per-peak substructure refinement...")
        refine_result = refine_peak_substructure_once(
            comps=comps,
            reads_df=reads_df,
            peak_subsets=peak_subsets,
            min_reads_for_refinement=min_reads_for_subpeak_refinement,
            min_child_fraction=min_subpeak_fraction,
            min_subpeak_distance=min_subpeak_distance,
            max_subpeak_sd=max_subpeak_sd,
            min_bic_gain_for_split=min_bic_gain_for_subpeak_split,
        )
        comps = refine_result.comps
        reads_df = refine_result.reads_df
        peak_subsets = refine_result.peak_subsets

    if comps.empty:
        finalize_no_itd("No GMM peaks remained after refinement.")

    total_effective_reads = int(comps.get("effective_read_count", pd.Series(dtype=int)).fillna(0).sum())
    if reads_df.empty or total_effective_reads == 0:
        finalize_no_itd("No reads remained assigned to a valid GMM peak after filtering.")

    logger.debug("Fitted GMM components (filtered):")
    logger.debug(comps)
    logger.debug("Per-read assignments (first 10 reads):")
    logger.debug(reads_df[:10])
    logger.debug("Subsets of reads per peak:")
    if logger.isEnabledFor(logging.DEBUG):
        for pid, reads in peak_subsets.items():
            logger.debug(f" Peak {pid} ({len(reads)} reads): {reads[:5]}{'...' if len(reads) > 5 else ''}")

    #Plot the size distribution and the identified peaks 
    plot_name = f"{sample_name}_itd_gmm_fit"
    out_prefix = os.path.join(flt3_data_folder, plot_name)
    plot_gmm_itds(
        reads_df=reads_df,
        comps=comps,
        bins=100,
        assign_mode=peak_read_assignment_mode,
        out_prefix=out_prefix,
        title="FLT3-ITD Read Length Distribution",
    )

    non_wt_comps = comps.loc[
        comps["peak_alias"].astype(str).str.upper() != "WT"
    ].copy()
    non_wt_comps = non_wt_comps.loc[non_wt_comps["effective_read_count"].fillna(0) > 0]
    if non_wt_comps.empty:
        finalize_no_itd("No ITD peaks were detected beyond the WT peak.")
        

    # Process per-peak reads and collect insertion calls.
    all_itd_insertions = []

    for alias, read_ids in peak_subsets.items():  # list of IDs per GMM peak
        if alias.upper() == "WT":
            continue
        if not read_ids:
            logger.info(f"Skipping ITD peak {alias}: no reads assigned after filtering.")
            continue
        logger.info(f"Processing ITD peak: {alias} ({len(read_ids)} reads)")
        df_itd = extract_itd_insertions_from_subset_parallel(
            reads_df=reads_df,                 # full read table
            read_ids_subset=read_ids,          # subset of read IDs for this ITD
            ref_seq=ref_seq,
            peak_alias=alias,
            comps=comps,
            threads=threads,
            itd_sd_factor=1.5
        )
        logger.debug(df_itd)
        all_itd_insertions.append(df_itd)

    if not all_itd_insertions:
        finalize_no_itd("No ITD insertions were detected after peak processing.")
    insertions_df = pd.concat(all_itd_insertions, ignore_index=True)
    logger.debug("Per-read insertion found (first 20 reads):")
    logger.debug(insertions_df[:20])

    # Save insertions table
    insertions_file = flt3_data_folder / f"{sample_name}_itd_insertions.tsv"
    insertions_df.to_csv(insertions_file, sep="\t", index=False)
    logger.info(f"Saved ITD insertions table: {insertions_file}")

    # Create ITD size distribution plot.
    logger.info("Creating ITD size distribution plot...")
    plot_itd_size_distribution(
        all_itd_insertions=insertions_df,
        out_dir=flt3_data_folder,
        sample_name=sample_name
    )

    # Build consensus sequences per peak using MSA
    logger.info("Building ITD consensus sequences per peak...")
    df_cons = build_itd_consensus_sequences(
    all_itd_insertions=insertions_df,
    comps=comps,
    out_dir=flt3_data_folder,
    sample_name=sample_name,
    max_unique=msa_max_unique,
    min_weight_coverage=msa_min_weight_coverage,
    base_threshold=msa_base_threshold,
    min_col_coverage=msa_min_col_coverage,
    ambiguous="N",
    )   

    logger.info("Per-peak consensus sequences:")
    logger.debug(df_cons)
    if df_cons.empty:
        finalize_no_itd("No ITD consensus sequences were generated.")

    # Build synthetic ITD references
    logger.info("Building ITD reference sequences and plots...")

    itd_refs = build_itd_reference_per_peak(
        df_cons=df_cons,
        ref_seq=ref_seq,
        genome_build=genome,
        amplicon_coords=amplicon_coords,
        exon_boundaries=exon_boundaries,
        out_dir=flt3_data_folder,
        sample_name=sample_name,
        exon_labels=exon_labels,
        flank_bp=200,
    )
    logger.info("Putative ITD reference sequences saved in the output folder")
    if not itd_refs:
        finalize_no_itd("No ITD reference sequences were generated.")

    # Create multi-FASTA for competitive alignment validation
    logger.info("Creating multi-reference FASTA for validation...")
    ref_dict, fasta_path = make_validation_refs(
    ref_seq=ref_seq,
    itd_ref_dict=itd_refs,
    sample_name=sample_name,
    out_dir=flt3_data_folder
    )

    # Prepare reads for validation
    logger.info("Preparing reads for validation alignment...")
    val_reads = prepare_validation_reads(
    reads_df=reads_df
    )

    if val_reads.empty:
        finalize_no_itd("No reads available for validation alignment.")

    # Align reads to multi-FASTA in parallel
    logger.info("Aligning reads to multi-reference FASTA in parallel...")
    validation_results = align_reads_multi_ref_parallel(val_reads, ref_dict, df_cons, logger=logger, threads=threads)
    if validation_results is None or validation_results.empty:
        finalize_no_itd("No validation alignments were produced.")
    logger.info("Validation alignment completed.")
    logger.debug("Validation alignment results (first 20 reads):")
    logger.debug(validation_results.iloc[0:20,3:])
    # Save validation results
    val_results_file = flt3_data_folder / f"{sample_name}_validation_read_support.tsv"
    validation_results.to_csv(val_results_file, sep="\t", index=False)
    logger.info(f"Saved validation results: {val_results_file}")

    # Calculate allele frequencies and strand bias
    logger.info("Calculating allele frequencies and strand bias...")
    summary_df = calculate_allele_frequencies_and_strand_bias(validation_results, min_allele_frequency=min_allele_frequency)
    if summary_df.empty:
        finalize_no_itd("No ITDs passed validation and allele-frequency filtering.")

    # Save VCF file with validated ITD calls
    logger.info("Writing validated ITD calls to VCF...")
    vcf_path = f"{sample_name}_FLT3_ITD_calls.vcf"
    output_path = os.path.join(output_folder, vcf_path)
    export_itd_vcf(summary_df, ref_dict, itd_refs, comps, output_path)

    if html_report:
        generate_itd_html_report(
            sample_name=sample_name,
            reference_genome=genome,
            seqio_reads=len(seqio_reads),
            val_reads=len(val_reads),
            summary_df=summary_df,
            itd_refs=itd_refs,
            output_dir=output_folder,
            plots_dir=flt3_data_folder,
            )

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
