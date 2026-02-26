# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

args = commandArgs(trailingOnly=TRUE)


# --- Load libraries -----------------------------------------------------------
library(dada2)

# --- Load data ----------------------------------------------------------------
datasets = c(args[1])
input.path = "TBD"
output.path = "TBD"

## Local function
ps.processing = function (dataset, path1 = input.path, path2 = input.path2, pattern.name = "_") {

  ## Read seqtab
  seqtab = readRDS(paste0(path2, dataset, "/", dataset, "_seqtab.rds"))
  
  ## Read taxa
  taxa = readRDS(paste0(path2, dataset, "/", dataset, "_taxa.rds"))
  
  ## Check order
  print(paste0("Samples order: ", all.equal(colnames(seqtab), rownames(taxa))))
  
  # Parse metadata
  samples = row.names(seqtab)
  metadata = read.delim(paste0(path1, dataset, "/00_parsed_metadata.tsv"))
  metadata$Dataset = dataset
  metadata$Run = metadata$Run
  rownames(metadata) = metadata$Run
  rownames(metadata) = metadata$Run 
  
  ## Parse seqtab
  rownames(seqtab) = str_split_fixed(rownames(seqtab), pattern = pattern.name, n = 2)[,1]
  
  ## Order seqtab and metadata (mantain seqtab order)
  metadata = metadata[which(metadata$Run %in% rownames(seqtab)),] 
  metadata = metadata[rownames(seqtab),]
  
  print(paste0("Metadata order: ", all.equal(rownames(metadata), rownames(seqtab))))
  
  ## Process NAs from taxa annotation
  taxa.temp = as.data.frame(taxa)
  taxa.temp = taxa.temp %>% select(-Species)
  
  for (i in c(1:nrow(taxa.temp))) {
    ## Select row
    x = taxa.temp[i,]
    
    if (is.na(x[,"Genus"]) == FALSE)  {
      annot = x
      
    } else if (is.na(x[,"Family"]) == FALSE) {
      annot = c(x[,!colnames(x) %in% c("Genus")], paste0(x[,"Family"], "_unclassified"))
      
      
    } else if (is.na(x[,"Order"]) == FALSE) {
      annot = c(x[,!colnames(x) %in% c("Genus", "Family")], rep(paste0(x[,"Order"], "_unclassified"), 2))
      
    } else if (is.na(x[,"Class"]) == FALSE) {
      annot = c(x[,!colnames(x) %in% c("Genus", "Family", "Order")], rep(paste0(x[,"Class"], "_unclassified"), 3))
      
    } else if (is.na(x[,"Phylum"])  == FALSE) {
      annot = c(x[,!colnames(x) %in% c("Genus", "Family", "Order", "Class")], rep(paste0(x[,"Phylum"], "_unclassified"), 4))
      
    } else if (is.na(x[,"Kingdom"])  == FALSE) {
      annot = c(x[,!colnames(x) %in% c("Genus", "Family", "Order", "Class", "Phylum")], rep(paste0(x[,"Kingdom"], "_unclassified"), 5))
    }
    
    ## Add row
    taxa.temp[i,] = unlist(annot)
    
  }
  
  print(paste0("Number of different genera: ", length(unique(taxa.temp[,"Genus"]))))
  taxa.temp = as.matrix(taxa.temp)
  
  ## Phyloseq object
  ps = phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), # axa_are_rows=FALSE porque están en las columnas
                sample_data(metadata), 
                tax_table(taxa.temp))
  
  ps.glom = tax_glom(ps, taxrank = "Genus")
  
  return(ps.glom) 
}


################################################################################
#   Taxonomic annotations                                                      #
################################################################################

for (dataset in datasets) {

  print(paste0("Executing ", dataset, " dataset"))

  ## Read data
  seqtab = readRDS(paste0(output.path, dataset, "/", dataset, "_seqtab.rds"))
  
  print(paste0("Start GTDB"))
  
  ## Assing GTBD taxonomy
  taxa = assignTaxonomy(seqtab, paste0(input.path, "GTBD/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz"), tryRC = TRUE, multithread=TRUE, verbose = TRUE)
  
  ## Write GTDB results
  write.table(taxa, paste0(output.path, dataset, "/", dataset, "_taxa.txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
  saveRDS(taxa, file = paste0(output.path, dataset, "/", dataset, "_taxa.rds"))
  
  print(paste0("End GTDB"))
}

################################################################################
#   Processing phyloseq object                                                 #
################################################################################
ps = ps.processing(dataset = "TBD", path1 = path1, path2 = path2)
