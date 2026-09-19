#!/usr/bin/env Rscript
# Denoise one peak's reads into ASVs with DADA2, and report which ASV each read
# was assigned to.
#
# Called by the 'dada2' backend of haplotype_split.py. Reads a FASTQ of a single
# GMM length peak and writes a two-column TSV (read_id, cluster_id) plus a FASTA
# of the inferred ASVs.
#
# DADA2's error model is built for Illumina substitutions. The settings below are
# the ones DADA2 documents for long, indel-prone reads (PacBio CCS), which is the
# closest supported regime to ONT:
#   BAND_SIZE = 32             wider alignment band, for indels
#   HOMOPOLYMER_GAP_PENALTY    softer gaps inside homopolymers, where ONT errors
#                              concentrate
#   USE_QUALS = TRUE           quality scores carry real signal here
#
# Usage:
#   Rscript dada2_cluster.R <in.fastq> <out_clusters.tsv> <out_asvs.fasta> \
#       [omega_a] [band_size] [homopolymer_gap_penalty] [min_asv_reads]

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: dada2_cluster.R <in.fastq> <out_clusters.tsv> <out_asvs.fasta> [omega_a] [band_size] [hp_gap_penalty] [min_asv_reads]")
}
fastq_in    <- args[1]
clusters_out<- args[2]
asvs_out    <- args[3]
omega_a     <- if (length(args) >= 4) as.numeric(args[4]) else 1e-40
band_size   <- if (length(args) >= 5) as.integer(args[5]) else 32L
hp_penalty  <- if (length(args) >= 6) as.numeric(args[6]) else -1
min_reads   <- if (length(args) >= 7) as.integer(args[7]) else 20L

setDadaOpt(BAND_SIZE = band_size, HOMOPOLYMER_GAP_PENALTY = hp_penalty)

# --- dereplicate -------------------------------------------------------------
derep <- derepFastq(fastq_in, verbose = FALSE)
n_reads <- sum(derep$uniques)
message(sprintf("[dada2_cluster] %d reads, %d unique", n_reads, length(derep$uniques)))

if (length(derep$uniques) < 2) {
  # nothing to separate
  writeLines(character(0), clusters_out)
  writeLines(character(0), asvs_out)
  quit(save = "no", status = 0)
}

# --- error model, learned from this peak -------------------------------------
# Learning on the peak itself rather than a global model: each peak is one
# amplicon at one length, so the error profile is homogeneous and the sample is
# what we actually want to denoise.
err <- tryCatch(
  learnErrors(fastq_in, multithread = TRUE, verbose = FALSE,
              errorEstimationFunction = loessErrfun, randomize = FALSE),
  error = function(e) {
    message("[dada2_cluster] learnErrors failed (", conditionMessage(e),
            "); falling back to a self-consistent estimate from dada()")
    NULL
  }
)

dd <- dada(derep, err = err, multithread = TRUE, verbose = FALSE,
           OMEGA_A = omega_a, selfConsist = is.null(err))

asv_seqs <- dd$sequence
message(sprintf("[dada2_cluster] %d ASVs inferred", length(asv_seqs)))

# --- map every read back to its ASV -----------------------------------------
# dd$map gives, for each unique sequence in derep, the index of the ASV it was
# assigned to. derep$map gives, for each input read, the index of its unique
# sequence. Composing the two gives read -> ASV.
read_to_unique <- derep$map
unique_to_asv  <- dd$map
read_to_asv    <- unique_to_asv[read_to_unique]

ids <- as.character(ShortRead::id(ShortRead::readFastq(fastq_in)))
ids <- sub("\\s.*$", "", ids)   # FASTQ id is the first whitespace-delimited token

keep <- !is.na(read_to_asv)
df <- data.frame(
  cluster_id = paste0("ASV", read_to_asv[keep]),
  read_id    = ids[keep],
  stringsAsFactors = FALSE
)

# drop ASVs below the minimum read count; the Python side folds these back
sizes <- table(df$cluster_id)
big <- names(sizes)[sizes >= min_reads]
message(sprintf("[dada2_cluster] %d/%d ASVs have >= %d reads",
                length(big), length(unique(df$cluster_id)), min_reads))

write.table(df, clusters_out, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

fa <- character(0)
for (i in seq_along(asv_seqs)) {
  n <- sum(df$cluster_id == paste0("ASV", i))
  fa <- c(fa, sprintf(">ASV%d size=%d", i, n), asv_seqs[i])
}
writeLines(fa, asvs_out)
