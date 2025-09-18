# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# PhD thesis: Sex differences and cross-disease molecular mechanisms in multiple 
# sclerosis: insights from transcriptomics and metagenomics analyses
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
set.seed(1234)

# --- Load libraries -----------------------------------------------------------
library(phyloseq)
library(DirichletMultinomial)
library(dplyr)
library(MetBrewer)
library(stringr)
library(GGally)

## Local functions
parsing2plot = function(genus_table, taxa_list, et_output, et_column, col.names = c("k1", "k2", "k3", "k4")){
  
  ## Change Genus matrix to subselect relevant taxa and merge all others in "Other"
  taxa.counts = data.frame(genus_table[, which(colnames(genus_table)%in%taxa_list)])
  other.taxa.counts = data.frame("Other" = rowSums(genus_table[,-which(colnames(genus_table) %in% taxa_list)]))
  
  taxa_table = data.frame(merge(taxa.counts, other.taxa.counts, by="row.names"), row.names=1)
  
  ## Merge Genus matrix with Enterotype definition by k4
  taxa_table = data.frame(merge(taxa_table, et_output[,et_column,drop=F], by="row.names"), row.names=1)
  taxa_table[,et_column] = as.factor(taxa_table[,et_column])
  
  ## Calculate mean expression
  plot_mat = t(data.frame(aggregate(. ~ get(et_column), taxa_table, mean), row.names=1))
  plot_mat = plot_mat[-nrow(plot_mat),] ## Eliminar el nombre de los clusters
  
  total = colSums(plot_mat)
  plot_mat_rel = t(plot_mat)
  plot_mat_rel = plot_mat_rel / total
  
  df = data.frame(t(plot_mat_rel))
  df$Bacteria = rownames(df)
  colnames(df) = c(col.names, "Bacteria")
  
  df2 = reshape::melt(df)
  
  return(df2)
  
}


plot.coord = function(ps, coord, variable, c.variable) {
  
  (p = plot_ordination(ps, coord, color = variable) +
     geom_point(size = 2) +
     stat_ellipse(type = "t", linetype = 3, geom = "polygon", aes(fill = variable),  alpha = 0.02) +
     theme_bw() +
     tune::coord_obs_pred() + # Same X and Y scale
     scale_color_manual(values = c.variable) +
     scale_fill_manual(values = c.variable) +
     theme(text = element_text(size = 16)))
  
  return(p)
}

# --- Load data ----------------------------------------------------------------

input.path = "TBD"
dataset = "TBD"
output.path = "TBD"

ps = readRDS(TBD)

################################################################################
#   Parsing data                                                               #
################################################################################
otu_table = t(data.matrix(otu_table(ps)))
otu_table = as.matrix(otu_table) 

otu_table = t(otu_table)

head(colnames(otu_table)) # Taxa
head(rownames(otu_table)) # Samples

dim(otu_table)

################################################################################
#   DMM Modelling                                                              #
################################################################################
all_dmns = 6 
dmn_list = numeric(all_dmns)

for(i in 1:all_dmns){
  print(i)
  assign(paste0("dmn_", i), dmn(otu_table, i, verbose=TRUE))
}

## List all DMN
dmn_list = list(dmn_1, dmn_2, dmn_3, dmn_4, dmn_5, dmn_6)

# Minimum number of Dirichlet Components
lplc = sapply(dmn_list, laplace)
BIC = sapply(dmn_list, BIC)
AIC = sapply(dmn_list, AIC)
dmn_list[[which.min(lplc)]] 
dmn_list[[which.min(BIC)]]
dmn_list[[which.min(AIC)]] 

# Plot
par(mfrow=c(1,3))
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit, Laplace")
plot(BIC, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit, BIC") #4 #optimal #minimal at plot = maximal information 
plot(AIC, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit, AIC") 


################################################################################
#   DMM MODELLING                                                              #
################################################################################
lplc = sapply(dmn_list, laplace)

dmn1 = dmn_list[[1]]
DM_1 = mixture(dmn1, assign = TRUE)

dmn2 = dmn_list[[2]]
DM_2 = mixture(dmn2, assign = TRUE)

dmn3 = dmn_list[[3]]
DM_3 = mixture(dmn3, assign = TRUE)

dmn4 = dmn_list[[4]]
DM_4 = mixture(dmn4, assign = TRUE)

dmn5 = dmn_list[[5]]
DM_5 = mixture(dmn5, assign = TRUE)

dmn6 = dmn_list[[6]]
DM_6 = mixture(dmn6, assign = TRUE)

DM_all = cbind(DM_1, DM_2, DM_3, DM_4, DM_5, DM_6)
colnames(DM_all) = c("k1","k2", "k3", "k4", "k5", "k6")

write.table(DM_all, paste0(output.path, dataset, "_results.tsv"), sep= "\t",col.names=NA)


################################################################################
#   Explore clusters                                                           #
################################################################################
et_column = "k4"

et_output = DM_all
table(et_output[,et_column])

samples = row.names(et_output)
et_output = data.frame(apply(et_output,2,function(x) as.factor(x)))
row.names(et_output) = samples

## Taxa to plot
genus_table=data.frame(otu_table) 
colnames(genus_table) = tax_table(ps@tax_table[,"Genus"])

### PRJNA684124
taxa_list=c("Bacteroides", "Bacteroides_F", 
            "Phocaeicola", "Phocaeicola_A",
            "Faecalibacterium",
            "Prevotella", 
            "Alistipes", "Alistipes_A",
            "Ruminococcus", "Ruminococcus_B", "Ruminococcus_C", "Ruminococcus_E", "Ruminococcus_F", "Ruminococcus_G",
            "Methanobrevibacter_A", "Methanobrevibacter_C",
            "Akkermansia") #TAXA to plot at end ("Bacteroides","Faecalibacterium","Prevotella" are mandatory)

colorss = c("Bacteroides" = "#ffd166", "Bacteroides_F" = "#ffd166", "Phocaeicola" = "#ffe6a7", "Phocaeicola_A" = "#ffe6a7", 
            "Faecalibacterium" = "#06d6a0",
            "Prevotella" = "#118ab2",
            "Alistipes" = "#fb5607", "Alistipes_A" = "#fb5607",
            "Ruminococcus" = "#a663cc", "Ruminococcus_B" = "#a663cc", "Ruminococcus_C" = "#a663cc", "Ruminococcus_E" = "#a663cc", "Ruminococcus_F" = "#a663cc", "Ruminococcus_G" = "#a663cc",
            "Methanobrevibacter_A" = "#a8dadc", "Methanobrevibacter_C" = "#a8dadc",
            "Akkermansia" = "#fbc4ab",    "Other" = "#e9ecef")

## Create data.frame with the clusters
ET_tab = genus_table[, et_column,drop=FALSE]
colnames(ET_tab) = "ET_clusters"

plot_mat = genus_table
plot_mat = t(data.frame(aggregate(. ~ k4,plot_mat,mean),row.names=1))

total = colSums(plot_mat)
plot_mat_rel = t(plot_mat)
plot_mat_rel = plot_mat_rel / total

df = data.frame(t(plot_mat_rel))
df$Bacteria = rownames(df)
colnames(df) = c("k1", "k2", "k3", "k4", "Bacteria")

df2 = reshape::melt(df)

df2$Bacteria = factor(df2$Bacteria, levels = rev(names(colorss)))
ggplot(df2, aes(x = variable, y = value, fill = Bacteria)) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(data=subset(df2, df2$value > 0.001), aes(label = paste0(substr(Bacteria, 1, 5), "   ", round(value, 2))),
            position = position_stack(vjust =0.5),
            size = 2.5,
            color = "black") +
  theme_bw()+ 
  ylab("Abundance ratio (from mean abudance)") + xlab("Cluster") +
  scale_fill_manual(values = colorss, name = "Selected bacteria")




