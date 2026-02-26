# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

args = commandArgs(trailingOnly=TRUE)

# --- Load libraries -----------------------------------------------------------
library(dada2)
library(dplyr)
library(stringr)
library(ggplot2)
library(ShortRead)

## Local functions
source("./01_dada2_functions.R")

# --- General variables --------------------------------------------------------

## Paths
setwd("TBD")
input.path = "TBD"
output.path = "TBD"

## Arguments
dataset = args[1] ## Dataset to evaluate
c.left = c(as.numeric(args[2]), as.numeric(args[3])) # Left trimming
c.right = c(as.numeric(args[4]), as.numeric(args[5])) # Right trimming
c.quim = as.numeric(args[6]) # Quimera value


################################################################################
#   Filtering and trimming                                                     #
################################################################################

results = order.seq.paired(input.path, dataset = dataset)

# Extract results
fq.forward = results[[1]]
fq.reverse = results[[2]]
samples = results[[3]]
rm(results)

# Plot individual quality profile of the interested read
pdf(paste0(output.path, dataset, "/quality_reads_dada2.pdf"), onefile=T) 
plotQualityProfile(fq.forward[1:2])
plotQualityProfile(fq.reverse[1:2])
dev.off()


# Name of filtered sequences
filt.forward = file.path(paste0(input.path, dataset, "/filtered"), paste0(samples, "_1_filt.fastq.gz"))
filt.reverse = file.path(paste0(input.path, dataset, "/filtered"), paste0(samples, "_2_filt.fastq.gz"))

# Paramters for trimming and filtering
results.trim = filter.trim(fqf = fq.forward, fqf_filtered = filt.forward,
                      fqr = fq.reverse, fqr_filtered = filt.reverse,
                      trimleft = c.left, trimright = c.right, minlen = 125, trunq = 2)

rm(filt.forward, filt.reverse)


# Load oredered data 
results = order.filt.seq.paired(input.path, dataset = dataset) # Ordered again just in case

# Extract results
filt.forward = results[[1]]
filt.reverse = results[[2]]
samples = results[[3]]
rm(results)

################################################################################
#   DADA2                                                                      #
################################################################################
results.denoised = all.denoising(filt.forward = filt.forward, filt.reverse = filt.reverse, output.path = output.path, dataset = dataset,
                  cutoff.quimera = c.quim)

saveRDS(results.denoised, file = paste0(output.path, dataset, "/results_denoised.rds"))

derep.forward = results.denoised[["derep.forward"]]
derep.reverse = results.denoised[["derep.reverse"]]

err.forward = results.denoised[["err.forward"]]
err.reverse = results.denoised[["err.reverse"]]

dada.forward = results.denoised[["dada.forward"]]
dada.reverse = results.denoised[["dada.reverse"]]

merge.reads = results.denoised[["merge.reads"]]
merge.nochim = results.denoised[["merge.nochim"]]

rm(results.denoised)

################################################################################
#   Statistics                                                                 #
################################################################################

# Track removed reads in each step
getN = function(x) sum(getUniques(x))

results.trim = results.trim[results.trim[,"reads.out"] != 0,]

track = as.data.frame(cbind(results.trim , sapply(dada.forward, getN), 
                            sapply(dada.reverse, getN), sapply(merge.reads, getN), 
                            sapply(merge.nochim, getN)))

colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track$ratio.filtered = track$filtered / track$input
track$ratio.merged =track$merged / track$filtered
track$ratio.nonchim = track$nonchim / track$merged

write.table(track, paste0(output.path, dataset, "/", dataset, "__stats.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)


################################################################################
#   ASV COUNT MATRIX                                                           #
################################################################################
seqtab = makeSequenceTable(merge.nochim)
dim(seqtab)

saveRDS(seqtab, file = paste0(output.path, dataset, "/", dataset, "_seqtab.rds"))



