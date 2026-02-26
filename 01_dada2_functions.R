# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #



# --- Local functions ----------------------------------------------------------

################################################################################
#   Order raw sequences                                                        #
################################################################################

order.seq.paired = function(input.path, dataset, Fpattern = "_1.fastq.gz", Rpattern = "_2.fastq.gz") {
  
  # Forward reads
  fqF = sort(list.files(paste0(input.path, dataset), pattern = Fpattern, full.names = TRUE))
  fqF = fqF[str_order(fqF, numeric = TRUE)]
  
  # Reverse reads
  fqR = sort(list.files(paste0(input.path, dataset), pattern = Rpattern, full.names = TRUE))
  fqR = fqR[str_order(fqR, numeric = TRUE)]
  
  # Sample names
  sample.names <- sapply(strsplit(basename(fqF), "_"), `[`, 1)
  
  return(list(fqF, fqR, sample.names))
}


################################################################################
#   Filtering and trimming                                                     #
################################################################################

filter.trim = function(fqf, fqf_filtered, fqr, fqr_filtered, trimleft = 0, trimright = 0, minlen, trunq, cutoff.trunc = 0){
  
  filtered = filterAndTrim(fqf, fqf_filtered, fqr, fqr_filtered, #truncLen = 275,
                           trimLeft = trimleft, trimRight = trimright, truncQ = trunq, minLen = minlen,
                           maxN=0, rm.phix = TRUE, truncLen = cutoff.trunc, #maxN lo dejamos a 0 porque dada2 no permite N
                           multithread = 1)
  return(filtered)
}


################################################################################
#   Order reads                                                                #
################################################################################

order.filt.seq.paired = function(input.path, dataset, Fpattern = "_1_filt.fastq.gz", Rpattern = "_2_filt.fastq.gz") {
  
  # Forward reads
  fqF <- sort(list.files(paste0(input.path, dataset, "/filtered"), pattern = Fpattern, full.names = TRUE))
  fqF <- fqF[str_order(fqF, numeric = TRUE)]
  
  # Reverse reads
  fqR <- sort(list.files(paste0(input.path, dataset, "/filtered"), pattern = Rpattern, full.names = TRUE))
  fqR <- fqR[str_order(fqR, numeric = TRUE)]
  
  # Sample names
  sample.names <- sapply(strsplit(basename(fqF), "_"), `[`, 1)
  
  return(list(fqF, fqR, sample.names))
}


################################################################################
#   Denoising                                                                  #
################################################################################

all.denoising = function(filt.forward, filt.reverse, output.path, dataset, cutoff.quimera) {
  
  # Dereplicate
  derep.forward = derepFastq(filt.forward, verbose=TRUE)
  derep.reverse = derepFastq(filt.reverse, verbose=TRUE)
  
  # Learn errors rate
  ## Forward reads
  err.forward = learnErrors(derep.forward, randomize=TRUE, multithread=FALSE)
  
  ### Plot F errors
  p = plotErrors(err.forward, nominalQ=TRUE) + ggtitle("Error frequencies for forward reads")
  
  svg(filename = paste0(output.path, dataset, "/", dataset, "_forward_errors.svg"), width = 10, height = 8)
  plot(p)
  dev.off()
  
  ## Reverse reads
  err.reverse = learnErrors(derep.reverse, randomize=TRUE, multithread=FALSE)
  
  ### Plot R errors
  p = plotErrors(err.reverse, nominalQ=TRUE) + ggtitle("Error frequencies for reverse reads")

  svg(filename = paste0(output.path, dataset, "/", dataset, "_reverse_errors.svg"), width = 10, height = 8)
  plot(p)
  dev.off()
  
  # Sample inference
  dada.forward = dada(derep = derep.forward, err = err.forward, multithread=FALSE) 
  dada.reverse = dada(derep = derep.reverse, err = err.reverse, multithread=FALSE)
  
  # Merge reads
  merge.reads = mergePairs(dada.forward, derep.forward, dada.reverse, derep.reverse, verbose = TRUE)
  
  # Remove quimeras (bimeras)
  merge.nochim = removeBimeraDenovo(merge.reads, multithread=FALSE, verbose=TRUE, method="consensus", minFoldParentOverAbundance=cutoff.quimera)
  
  # All results
  all.results = list(derep.forward, derep.reverse, err.forward, err.reverse, dada.forward, dada.reverse, merge.reads, merge.nochim)
  names(all.results) = c("derep.forward", "derep.reverse", "err.forward", "err.reverse", "dada.forward", "dada.reverse", "merge.reads", "merge.nochim")
  
  # Return results
  return(all.results)
  
}


