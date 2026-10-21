suppressPackageStartupMessages(library(dada2))
p <- '/workspace/review/fastq_51_57_analysis'
fq <- file.path(p,'haplotypes/dada2_ITD_1/reads.fastq')
setDadaOpt(BAND_SIZE=32,HOMOPOLYMER_GAP_PENALTY=-1)
dr <- derepFastq(fq)
err <- learnErrors(fq,multithread=4,verbose=FALSE,errorEstimationFunction=loessErrfun,randomize=FALSE)
for (omega in c(1e-40,1e-10,1e-3)) {
 dd <- dada(dr,err=err,multithread=4,verbose=FALSE,OMEGA_A=omega)
 cat('OMEGA',omega,'ASVS',length(dd$sequence),'LENGTHS',nchar(dd$sequence),'COUNTS',dd$abundance,'\n')
 write.table(dd$clustering,file.path(p,paste0('dada2_clustering_',omega,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 if (omega==1e-40) {
  cat('PACKAGE_VERSION',as.character(packageVersion('dada2')),'\n')
  # Diagnostic p-values for the frequent, exact 393 bp putative ITD2 read.
  s <- names(dr$uniques)[nchar(names(dr$uniques))==393]
  target <- s[which.max(dr$uniques[s])]
  cat('ITD2_UNIQUE_COUNT',dr$uniques[target],'\n')
  print(dd$birth_subs)
  saveRDS(dd,file.path(p,'dada2_default.rds'))
 }
}
