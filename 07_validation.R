# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Load libraries --------------------------------------------------------- #
library(phyloseq)
library(DESeq2)
library(xlsx)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(PCAtools)
library(pheatmap)
library(microbiome)
library(MetBrewer)
library(dplyr)
library(DT)
library(rstatix)
library(boot)
library(rstatix)
library(boot)

# --- Local functions ---------------------------------------------------------

taxa.filtering.func = function(ps, prev.threshold = prev.value, sample.threshold = sample.value) {
  
  ## Relative counts
  ps.rel = transform_sample_counts(ps, function(x) x / sum(x))
  
  ## Extract otu table
  otu.tab = as.data.frame(otu_table(ps.rel))
  if (taxa_are_rows(ps.rel)) {
    otu.tab = t(otu.tab)
  }
  colnames(otu.tab) = ps@tax_table[,"Genus"]
  
  ## Check thresholds
  keep_taxa = apply(otu.tab, 2, function(taxon_abund) {
    (sum(taxon_abund > prev.threshold) / length(taxon_abund)) >= sample.threshold
  })
  
  table(keep_taxa)
  
  return(colnames(otu.tab)[keep_taxa])
}


# --- Load data ---------------------------------------------------------------
input.path = "TBD"
PRJEB32762 = readRDS("TBD")

dim(PRJEB32762@otu_table)
taxa_names(PRJEB32762) = tax_table(PRJEB32762)[,"Genus"]

## Variable to control
PRJEB32762@sam_data$to_correct = paste0(PRJEB32762@sam_data$Cohort.y, "_", PRJEB32762@sam_data$site)
table(PRJEB32762@sam_data$to_correct)

ps = PRJEB32762


################################################################################
#   Filter low abundant taxa                                                   #
################################################################################
taxa.dataset = taxa.filtering.func(ps = PRJEB32762, prev.threshold = TBD, sample.threshold = TBD) 
PRJEB32762.filt = prune_taxa(taxa.dataset, PRJEB32762) 

################################################################################
#   Normalization                                                              #
################################################################################
counts = t(as.data.frame(PRJEB32762.filt@otu_table))

## Normalized
PRJEB32762.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))
PRJEB32762.vst = PRJEB32762.vst[-nrow(PRJEB32762.vst),]

################################################################################
#   Diferential abundance analysis                                             #
################################################################################

## Prepare data
dataset.normalized = t(PRJEB32762.vst)
dataset.metadata = as.data.frame(sample_data(PRJEB32762.filt))
dataset.name = "PRJEB32762"

## Parsing data
df = merge(x = dataset.normalized, by.x = "row.names", y = dataset.metadata, by.y = "row.names")
rownames(df) = df$Row.names
df = df[,-1]
all.tests = gsub(pattern = "-", replacement = "", taxa.dataset)

## Convert as factors
df$Condition = factor(df$Condition, levels = c("Control", "MS"))
df$Sex = factor(df$Sex, levels = c("Male", "Female"))

comparisons = c("Control", "Female", "Male", "MS")
path = ""

for (comparison in comparisons) {
  
  print(comparison)
  significant.MA = readRDS(paste0("TBD", folder.name, "/", comparison, "/significant_results_MA.rds"))
  
  if(TRUE == TRUE) { #nrow(significant.MA > -1)
    
    print(paste0("Significant results in MA: ", nrow(significant.MA)))
    print(paste0((significant.MA$Taxa)))
    print(paste0((significant.MA$Taxa), collapse = ", "))
    
    ## Taxa to validate
    taxa2validate = significant.MA$Taxa
    print(table(taxa2validate %in% all.tests))
    tests = taxa2validate[taxa2validate %in% all.tests]

    if (length(tests) > 0) {
    
      # ----------------------- #
      # ---- Wilcoxon test ---- #
      # ----------------------- #
      
      if (comparison %in% c("Control", "MS")) {
        
        df.test = df %>% filter(Condition == comparison)
        results.w = wilcoxon.latent.local(data = df.test, vars.test = tests, var.group = "Sex", var.latent = "to_correct")
        
        results.w$Comparison = comparison
        
        df.temp = data.frame(matrix(ncol = 6, nrow = 0))
        colnames(df.temp) = c("mean1", "mean2", "effsize", "n1", "n2", "magnitude")
        
        for (z in 1:nrow(results.w)) {
        
          genus = results.w$Variable[z]
          
          mean1 = mean(df.test[df.test$Sex == "Male", genus]) ## Group 1 reference
          mean2 = mean(df.test[df.test$Sex == "Female", genus])
          
          formula = as.formula(paste0(genus, " ~ Sex"))
          
          effect.size =  wilcox_effsize(data = df.test,
                                        formula = formula, comparisons = list(c("Female", "Male")), ci = FALSE) 
          
          ## Determine sign of effect size
          if (mean1 > mean2) { # Group 1 reference
            effect.size$effsize = effect.size$effsize * -1
          }
          df.temp = rbind(df.temp, c("mean1" = mean1, "mean2" = mean2, effect.size[4:7]))
        }
        
        all.results = cbind(results.w, df.temp)
        
      }
      
      if (comparison %in% c("Female", "Male")) {
        
        df.test = df %>% filter(Sex == comparison)
        results.w = wilcoxon.latent.local(data = df.test, vars.test = tests, var.group = "Condition", var.latent = "to_correct")
        
        results.w$Comparison = comparison
        
        df.temp = data.frame(matrix(ncol = 6, nrow = 0))
        colnames(df.temp) = c("mean1", "mean2", "effsize", "n1", "n2", "magnitude")
        
        for (z in 1:nrow(results.w)) {
          
          genus = results.w$Variable[z]
          
          mean1 = mean(df.test[df.test$Condition == "Control", genus]) ## Group 1 reference
          mean2 = mean(df.test[df.test$Condition == "MS", genus])
          
          formula = as.formula(paste0(genus, " ~ Condition"))
          
          effect.size =  wilcox_effsize(data = df.test,
                                        formula = formula, comparisons = list(c("MS", "Control")), ci = FALSE)    
          
          ## Determine sign of effect size
          if (mean1 > mean2) { # Group 1 reference
            effect.size$effsize = effect.size$effsize * -1
          }
          
          df.temp = rbind(df.temp, c("mean1" = mean1, "mean2" = mean2, effect.size[4:7]))
        }

        
      }
      
      all.results = cbind(results.w, df.temp)
      
      # ## Save results
      write.table(x = all.results, file = paste0(path, comparison, "_", "wilcoxon.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
      saveRDS(object = all.results, file = paste0(path, comparison, "_", "wilcoxon.rds"))

      print("Raw p.value")
      print(paste0("logFC > 0: ", nrow(all.results[all.results$p.value < 0.05 & all.results$effect.size > 0,])))
      print(paste0("logFC < 0: ", nrow(all.results[all.results$p.value < 0.05 & all.results$effect.size < 0,])))
      
      print("Adjusted p.value")
      print(paste0("logFC > 0: ", nrow(all.results[all.results$p.adjusted < 0.05 & all.results$effect.size > 0,])))
      print(paste0("logFC < 0: ", nrow(all.results[all.results$p.adjusted < 0.05 & all.results$effect.size < 0,])))
      
    }
  }
}



