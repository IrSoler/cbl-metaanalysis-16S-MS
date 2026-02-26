# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# --- Load libraries -----------------------------------------------------------
library(phyloseq)
library(ggplot2)

## Local functions
plot.taxa.local = function(ps) {
  
  taxa = data.frame(otu_table(ps))
  df = data.frame(sample_data(ps))
  
  taxa$Ntaxa = Matrix::rowSums(taxa > 0)
  taxa$Ncounts = Matrix::rowSums(taxa[,-length(taxa)])
  taxa$Sample = rownames(taxa)
  
  df2plot = merge(taxa, by.x = "row.names", df, by.y = "row.names")
  rownames(df2plot) = df2plot$Row.names
  
  colorss = c("TBD")

  # Plot N taxa
  ggplot() + geom_col(data = df2plot, mapping = aes(x = Sample, y = Ntaxa, fill = TBD)) + theme_bw() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 6)) +
    scale_fill_manual(values = colorss)
}


plot.counts.local = function(ps) {
  
  taxa = data.frame(otu_table(ps))
  df = data.frame(sample_data(ps))
  
  taxa$Ntaxa = Matrix::rowSums(taxa > 0)
  taxa$Ncounts = Matrix::rowSums(taxa[,-length(taxa)])
  taxa$Sample = rownames(taxa)
  
  df2plot = merge(taxa, by.x = "row.names", df, by.y = "row.names")
  
  colorss = c("TBD")

  # Plot N counts
  ggplot() + geom_col(data = df2plot, mapping = aes(x = Sample, y = Ncounts, fill = TBD)) + theme_bw() +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, size = 6)) +
    scale_fill_manual(values = colorss)
}


plot.rarefraction.local = function(ps) {
  data2rarefy = as(otu_table(ps), "matrix")
  class(data2rarefy)
  rarecurve(data2rarefy[40:47,], step = 1000, label = FALSE, cex = 0.5)
}


# --- Load data ----------------------------------------------------------------
ps.dataset = readRDS("TBD")

################################################################################
#   Exploratory plots                                                          #
################################################################################
plot.taxa.local(ps = ps.dataset)
plot.counts.local(ps = ps.dataset)
plot.rarefraction.local(ps = ps.dataset)
