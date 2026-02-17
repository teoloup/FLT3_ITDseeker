#!/usr/bin/env python3
"""
BAM Read Extractor Module
Extracts reads mapping to FLT3 region and performs primer trimming
"""

import logging
import pysam
from pathlib import Path
from typing import List, Dict, Tuple
from dataclasses import dataclass
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import tempfile
import subprocess
import os
import shutil


logger = logging.getLogger(__name__)


class FLT3ReadExtractor:
    """Extract and process reads from BAM file"""

    def __init__(self, bam_file: str, genome_build ,  min_mapping_quality, threads, min_length, max_itd_length, wt_amplicon_length, temp_dir, config=None ):
        self.bam_file = bam_file
        self.genome_build = genome_build
        self.min_mapping_quality = min_mapping_quality
        self.threads = threads
        self.min_length = min_length
        self.max_itd_length = max_itd_length
        self.temp_dir = temp_dir
        self.wt_amplicon_length = wt_amplicon_length

        if config is None:

            # Minimal fallback config with required attributes
            class MinimalConfig:
                forward_primer = "AGCAATTTAGGTATGAAAGCCAGC"
                reverse_primer = "CTGTACCTTTCAGCATTTTGACG" 
                start_hg38 = 28033301
                end_hg38 = 28034800
                start_hg19 = 28607438
                end_hg19 = 28608937
            self.config = MinimalConfig()
        else:
            self.config = config
        self.setup_coordinates()
        
    def setup_coordinates(self):
        """Set up genomic coordinates based on genome build"""
        if self.genome_build == "hg38":
            self.start = self.config.start_hg38
            self.end = self.config.end_hg38
        else:  # hg19
            self.start = self.config.start_hg19
            self.end = self.config.end_hg19
            
        # Check chromosome naming in BAM
        with pysam.AlignmentFile(self.bam_file, "rb") as bam:
            references = bam.references
            if "chr13" in references:
                self.chromosome = "chr13"
            elif "13" in references:
                self.chromosome = "13"
            else:
                raise ValueError("Cannot find chromosome 13 in BAM file")
                
        self.region = f"{self.chromosome}:{self.start}-{self.end}"
        logger.info(f"Using region: {self.region}")
    
    def extract_reads(self):
        """Extract reads mapping to FLT3 region using samtools view + fastq, the fastq file"""
        try:
            temp_dir = Path(str(self.temp_dir))  # Ensure temp_dir is a Path object
            logger.debug(f"Creating output folder: {temp_dir}")
            temp_dir.mkdir(parents=True, exist_ok=True)
            logger.debug(f"Output folder: {temp_dir}")
        except Exception as e:
            logger.error(f"Error creating temp folder: {e}")
            

        fastq_out = os.path.join(temp_dir, "region_reads.fastq")
        threads = self.threads
        region = f"{self.chromosome}:{self.start}-{self.end}"
        # Compose samtools view (region + quality filter) piped to samtools fastq
        samtools_view_cmd = [
            "samtools", "view",
            "-h",  # include header so samtools fastq works
            "--threads", str(threads), #threads to use
            "-F", "4",  # skip unmapped, secondary, supplementary
            "-q", str(self.min_mapping_quality),
            self.bam_file,
            region
        ]
        # Use shell redirection for samtools fastq output
        fastq_cmd_str = f"samtools fastq - > '{fastq_out}'"
        logger.info(f"Running samtools view + fastq for region: {' '.join(samtools_view_cmd)} | {fastq_cmd_str}")
        view_proc = subprocess.Popen(samtools_view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        fastq_proc = subprocess.Popen(fastq_cmd_str, shell=True, stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, fastq_stderr = fastq_proc.communicate()
        view_proc.stdout.close()
        view_stderr = view_proc.stderr.read()
        view_proc.stderr.close()
        view_proc.wait()
        if fastq_proc.returncode != 0 or view_proc.returncode != 0:
            logger.error(f"Samtools view/fastq failed: {view_stderr.decode()} {fastq_stderr.decode()}")
            shutil.rmtree(temp_dir)
            raise RuntimeError("Samtools view/fastq failed for region extraction.")
        logger.info(f"Extracted reads written to {fastq_out}")

        return fastq_out

    def trim_primers(self, fastq: str):
        """Trim primers using cutadapt with linked adapters and Nanopore error tolerance. Always keep trimmed FASTQ file for validation."""

        error_rate = 0.2
        threads = self.threads
        amplicon_length = self.wt_amplicon_length
        min_length = self.min_length
        fwd = self.config.forward_primer
        rev = str(Seq(self.config.reverse_primer).reverse_complement())

        logger.debug(f"Using {threads} threads for cutadapt")

        # Prepare temporary files
        #temp_dir = tempfile.mkdtemp(prefix="flt3_cutadapt_")
        fastq_in = fastq
        fastq_out = Path(self.temp_dir) / "trimmed.fastq"
        logger.debug(f"Trimmed FASTQ output: {fastq_out}")
        trimmed_fastq = str(fastq_out)

        # Prepare cutadapt command
        # Using linked adapters to find and trim everything OUTSIDE the amplicon
        # This keeps the amplicon (including primers) and removes flanking sequences
        cutadapt_cmd = [
            "cutadapt",
            "-g", f"{fwd}...{rev};rightmost",  # linked adapters: find forward then reverse
            "-e", str(error_rate),
            "-j", str(threads),  # number of parallel jobs
            "-m", str(min_length),  # min length after trimming
            "--rc",  # reverse complement search also
            "--discard-untrimmed",  # discard reads without both primers
            "--action=retain",  # keep primer sequences
            "--maximum-length", str(amplicon_length + int(self.max_itd_length * 1.1)),  # amplicon + max ITD + 10% buffer size
            "-o", trimmed_fastq,
            fastq_in
        ]

        logger.info(f"Running cutadapt for primer trimming: {' '.join(cutadapt_cmd)}")
        logger.debug(f"Full cutadapt command: {cutadapt_cmd}")  # Debug the actual command list
        result = subprocess.run(cutadapt_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logger.error(f"Cutadapt failed: {result.stderr.decode()}")
            shutil.rmtree(temp_dir)
            raise RuntimeError("Cutadapt failed for primer trimming.")

        # Parse trimmed FASTQ
        trimmed_records = []
        untrimmed_reads = list(SeqIO.parse(fastq_in, "fastq"))
        trimmed_reads = list(SeqIO.parse(trimmed_fastq, "fastq"))
        for rec in SeqIO.parse(trimmed_fastq, "fastq"):
            if rec.description.endswith(" rc"):
             rec.id = rec.id.replace(" rc", "")
             rec.seq = rec.seq.reverse_complement()
            trimmed_records.append(rec)

        reads = {rec.id: str(rec.seq) for rec in trimmed_records}

        logger.info(f"Cutadapt trimmed {len(trimmed_reads)} reads (input: {len(untrimmed_reads)})")
        # Do not delete temp_dir here; FASTQ is needed for validation
        return trimmed_fastq, reads
    
    def process_reads(self):
        """Extract reads and trim primers using cutadapt. Always keep trimmed FASTQ file for validation."""
        flt3_fastq = self.extract_reads()
        trimmed_flt3_fastq, trimmed_reads = self.trim_primers(flt3_fastq)

        return trimmed_reads, trimmed_flt3_fastq

def extract_flt3_reads(bam_file: str, genome_build: str , threads: int ,
                      min_mapping_quality: int , min_length: int , max_itd_length: int , wt_amplicon_length: int , temp_dir: str , config=None) -> Tuple[List[Dict], str]:
    """Main function to extract and process FLT3 reads, using main config if provided.
    Returns processed_reads and the trimmed FASTQ file path for centralized cleanup."""


    extractor = FLT3ReadExtractor(bam_file, genome_build, min_mapping_quality, threads = threads, min_length=min_length, max_itd_length=max_itd_length, wt_amplicon_length=wt_amplicon_length, temp_dir=temp_dir, config=config)
    processed_reads, fastq_file = extractor.process_reads()
    #shutil.rmtree(temp_dir) wait to decide wether i return a dict or the fastq, if debug keep temp
    return processed_reads, fastq_file