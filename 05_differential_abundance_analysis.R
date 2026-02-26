# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
set.seed(1234)

args = commandArgs(trailingOnly=TRUE)

# --- Arguments ----------------------------------------------------------------
prev.value = as.numeric(args[1])
sample.value = as.numeric(args[2])
folder.name = as.character(args[3])


# --- Load libraries -----------------------------------------------------------
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

## Local functions

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


wilcox_eff_stat_sex = function(data, indices, df.test2 = df.test, y = z) {
  
  # Subsampling
  d = data[indices, , drop = FALSE]  
  
  # Effect size for simulation
  res = wilcox_effsize(
    data = d,
    formula = as.formula(paste0(colnames(df.test2)[y], " ~ Sex"))
  )
  
  return(res$effsize) 
}


wilcox_eff_stat_condition = function(data, indices, df.test2 = df.test, y = z) {
  
  # Subsampling
  d = data[indices, , drop = FALSE]  
  
  # Effect size for simulation
  res = wilcox_effsize(
    data = d,
    formula = as.formula(paste0(colnames(df.test2)[y], " ~ Condition"))
  )
  
  return(res$effsize)
}


## Complete function
effect.size.func = function(df.func, group1.func, group2.func, Ntaxa.func) {
  
  df.effect.size = data.frame(matrix(ncol = 12, nrow = 0))
  
  for (z in 1:Ntaxa.func) {
    
    print(paste0("Processing: ", colnames(df.test)[z]))
    
    ## Calculate effect size
    if (comparison %in% c("Control", "MS")) {
      
      mean1 = mean(df.func[df.func$Sex == group1.func, z]) ## Group 1 reference
      mean2 = mean(df.func[df.func$Sex == group2.func, z])
      
      formula = as.formula(paste0(colnames(df.func)[z], " ~ Sex"))
      
      effect.size =  wilcox_effsize(data = df.func,
                                    formula = formula, comparisons = list(c("Female", "Male")), ci = FALSE)     
      
      boot_res = boot(
        data = df.func,
        statistic = wilcox_eff_stat_sex,
        y = z,
        df.test2 = df.func,
        R = 1000,
        strata = df.func$Sex
      )
      
      print(paste0("Number of NAs: ", sum(is.na(boot_res$t))))
      
      if(sum(is.na(boot_res$t)) == 1000) {
        
        effect.size$effsize = 0
        effect.size$magnitude = "small"
        effect.size$conf.low = 0
        effect.size$conf.high = 0
        
      } else {
      
      boot_ci_res = boot.ci(
        boot_res,
        conf = 0.95,
        type = "perc"
      )
      
      ci_low  = boot_ci_res$percent[ ,4]  # límite inferior
      ci_high = boot_ci_res$percent[ ,5]  # límite superior
      
      ## Add results to the data.frame
      effect.size$conf.low = ci_low
      effect.size$conf.high = ci_high
      }
    }
    
    if (comparison %in% c("Female", "Male")) {
      
      mean1 = mean(df.func[df.func$Condition == group1.func, z]) ## Group 1 reference
      mean2 = mean(df.func[df.func$Condition == group2.func, z])
      
      effect.size =  wilcox_effsize(data = df.func,
                                    formula = as.formula(paste0(colnames(df.func)[z], " ~ Condition")), comparisons = list(c("MS", "Control")), ci = FALSE)     
      
      boot_res = boot(
        data = df.func,
        statistic = wilcox_eff_stat_condition,
        y = z,
        df.test2 = df.func,
        R = 1000,
        strata = df.func$Condition
      )
      
      print(paste0("Number of NAs: ", sum(is.na(boot_res$t))))
      
      if(sum(is.na(boot_res$t)) == 1000) {
        
        effect.size$effsize = 0
        effect.size$magnitude = "small"
        effect.size$conf.low = 0
        effect.size$conf.high = 0
        
      } else {
      
      boot_ci_res = boot.ci(
        boot_res,
        conf = 0.95,
        type = "perc"
      )
      
      ci_low  = boot_ci_res$percent[ ,4]  # límite inferior
      ci_high = boot_ci_res$percent[ ,5]  # límite superior
      
      ## Add results to the data.frame
      effect.size$conf.low = ci_low
      effect.size$conf.high = ci_high
      
      }
    }
    
    effect.size$SE = (effect.size$conf.high - effect.size$conf.low) / (2 * 1.96)
    
    ## Determine sign of effect size
    if (mean1 > mean2) { # Group 1 reference
      effect.size$effsize = effect.size$effsize * -1
    }
    
    ## Parsing to provide results
    effect.size = as.data.frame(effect.size)
    effect.size$mean.group1 = mean1
    effect.size$mean.group2 = mean2
    
    print(paste0("Processed: ", colnames(df.test)[z]))
    
    df.effect.size = rbind(df.effect.size, effect.size)
    
  }
  
  return(df.effect.size)
  
}



# --- Load ps datasets -------------------------------------------------------- #
input.path = "TBD"

## Individual datasets
PRJNA684124 = readRDS(paste0(input.path, "PRJNA684124", "/", "PRJNA684124", "_ps.rds"))
dim(PRJNA684124@otu_table)
taxa_names(PRJNA684124) = tax_table(PRJNA684124)[,"Genus"]

PRJEB99111 = readRDS(paste0(input.path, "PRJEB99111", "/", "PRJEB99111", "_ps.rds"))
dim(PRJEB99111@otu_table)
taxa_names(PRJEB99111) = tax_table(PRJEB99111)[,"Genus"]

PRJNA732670 = readRDS(paste0(input.path, "PRJNA732670", "/", "PRJNA732670", "_ps.rds"))
dim(PRJNA732670@otu_table)
taxa_names(PRJNA732670) = make.unique(tax_table(PRJNA732670)[,"Genus"]) ## For UBA1381_unclassified (x2)

PRJEB34168 = readRDS(paste0(input.path, "PRJEB34168", "/", "PRJEB34168", "_ps.rds"))
dim(PRJEB34168@otu_table)
taxa_names(PRJEB34168) = tax_table(PRJEB34168)[,"Genus"]

PRJNA565173 = readRDS(paste0(input.path, "PRJNA565173", "/", "PRJNA565173", "_ps.rds"))
dim(PRJNA565173@otu_table)
taxa_names(PRJNA565173) = make.unique(tax_table(PRJNA565173)[,"Genus"]) ## For UBA1381_unclassified (x2)

PRJEB67783 = readRDS(paste0(input.path, "PRJEB67783", "/", "PRJEB67783", "_ps.rds"))
dim(PRJEB67783@otu_table)
taxa_names(PRJEB67783) = tax_table(PRJEB67783)[,"Genus"]


################################################################################
#   Filter low prevalent taxa                                                  #
################################################################################
ps.list = list(PRJNA684124, PRJEB99111, PRJNA732670, PRJEB34168, PRJNA565173, PRJEB67783)

taxa.dataset = lapply(ps.list, taxa.filtering.func)
common.taxa = Reduce(intersect, taxa.dataset)
print(paste0("N common taxa: ", length(common.taxa)))

## Filtered datasets
PRJNA684124.filt = prune_taxa(common.taxa, PRJNA684124) 
PRJEB99111.filt = prune_taxa(common.taxa, PRJEB99111) 
PRJNA732670.filt = prune_taxa(common.taxa, PRJNA732670) 
PRJEB34168.filt = prune_taxa(common.taxa, PRJEB34168) 
PRJNA565173.filt = prune_taxa(common.taxa, PRJNA565173) 
PRJEB67783.filt = prune_taxa(common.taxa, PRJEB67783) 

################################################################################
#   Normalization                                                              #
################################################################################

## Filtered datasets
PRJNA684124.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJNA684124.filt@otu_table))))
PRJEB99111.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJEB99111.filt@otu_table))))
PRJNA732670.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJNA732670.filt@otu_table))))
PRJEB34168.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJEB34168.filt@otu_table))))
PRJNA565173.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJNA565173.filt@otu_table))))
PRJEB67783.vst = DESeq2::varianceStabilizingTransformation(as.matrix(as.data.frame(t(PRJEB67783.filt@otu_table))))


################################################################################
#   Diferential abundance analyses                                             #
################################################################################

datasets.names = c("PRJNA684124", "PRJEB99111", "PRJNA732670", "PRJEB34168", "PRJNA565173", "PRJEB67783") 
datasets.filtered = list(PRJNA684124.filt, PRJEB99111.filt, PRJNA732670.filt, PRJEB34168.filt, PRJNA565173.filt, PRJEB67783.filt) ## PRJEB34168.filt
datasets.normalized = list(PRJNA684124.vst, PRJEB99111.vst, PRJNA732670.vst, PRJEB34168.vst, PRJNA565173.vst, PRJEB67783.vst) ## PRJEB34168.vst

Ntaxa = length(common.taxa)
comparisons = c("Control", "Female", "Male", "MS")
path = paste0("TBD", folder.name, "/")

for (i in 1:length(datasets.names)) {
  
  ## Extract data
  dataset.normalized = t(datasets.normalized[[i]])
  dataset.metadata = as.data.frame(sample_data(datasets.filtered[[i]]))
  dataset.name = datasets.names[i]
  
  ## Process names to avoid problems with the formula for the effect size
  colnames(dataset.normalized) = gsub(pattern = "-", replacement = "", colnames(dataset.normalized))
  
  ## Parsing data
  df = merge(x = dataset.normalized, by.x = "row.names", y = dataset.metadata, by.y = "row.names")
  rownames(df) = df$Row.names
  df = df[,-1]
  tests = gsub(pattern = "-", replacement = "", common.taxa)
  
  ## Convert as factors
  df$Condition = factor(df$Condition, levels = c("Control", "MS"))
  df$Sex = factor(df$Sex, levels = c("Male", "Female"))
  
  for (comparison in comparisons) {
    
    if (comparison == "Control") {
      
      ## Wilcoxon test
      df.test = df %>% filter(Condition == comparison)
      results.w = wilcoxon.local(data = df.test, vars.test = tests, var.group = "Sex")
      results.w.eff = effect.size.func(df.func = df.test, group1.func = "Male", group2.func = "Female", Ntaxa.func = Ntaxa) # Group 1 reference
      
      
      ## Define groups
      ordered.groups = df.test$Sex
      
      ## Define groups
      group1 = which(ordered.groups == "Male")
      group2 = which(ordered.groups == "Female")
      
      
    } else if(comparison == "Female") {
      
      ## Wilcoxon test
      df.test = df %>% filter(Sex == comparison)
      results.w = wilcoxon.local(data = df.test, vars.test = tests, var.group = "Condition")
      results.w.eff = effect.size.func(df.func = df.test, group1.func = "Control", group2.func = "MS", Ntaxa.func = Ntaxa) # Group 1 reference
      
      
      ## Define groups
      ordered.groups = df.test$Condition
      
      ## Define groups
      group1 = which(ordered.groups == "Control")
      group2 = which(ordered.groups == "MS")
      
      
    } else if(comparison == "Male") {
      
      ## Wilcoxon test
      df.test = df %>% filter(Sex == comparison)
      results.w = wilcoxon.local(data = df.test, vars.test = tests, var.group = "Condition")
      results.w.eff = effect.size.func(df.func = df.test, group1.func = "Control", group2.func = "MS", Ntaxa.func = Ntaxa) # Group 1 reference
      
      
      ## Define groups
      ordered.groups = df.test$Condition
      
      ## Define groups
      group1 = which(ordered.groups == "Control")
      group2 = which(ordered.groups == "MS")
      
      
    } else if (comparison == "MS") {
      
      ## Wilcoxon test
      df.test = df %>% filter(Condition == comparison)
      results.w = wilcoxon.local(data = df.test, vars.test = tests, var.group = "Sex")
      results.w.eff = effect.size.func(df.func = df.test, group1.func = "Male", group2.func = "Female", Ntaxa.func = Ntaxa) # Group 1 reference
      
      
      ## Define groups
      ordered.groups = df.test$Sex
      
      ## Define groups
      group1 = which(ordered.groups == "Male")
      group2 = which(ordered.groups == "Female")
      
    }
    
    
    ## Join with Wilcoxon result
    all.results = merge(x = results.w, by.x = "Variable", y = results.w.eff, by.y = ".y.")
    
    all.results = all.results[,-c(7,8)]
    
    colnames(all.results) = c("Taxa", "Group1", "Group2", "W", "p.value", "p.adjusted", "effect.size", "N1", "N2",  "magnitude", "CI.low", "CI.high", "SE", "mean.Group1", "mean.Group2")
    
    ## Save results
    write.table(x = all.results, file = paste0(path, dataset.name, "/", dataset.name, "_", comparison, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
    saveRDS(object = all.results, file = paste0(path, dataset.name, "/", dataset.name, "_", comparison, ".rds"))
    
    ## Print results
    print(dataset.name)
    print(comparison)
    
    print("Raw p.value")
    print(paste0("logFC > 0: ", nrow(all.results[all.results$p.value < 0.05 & all.results$effect.size > 0,])))
    print(paste0("logFC < 0: ", nrow(all.results[all.results$p.value < 0.05 & all.results$effect.size < 0,])))
    
    print("Adjusted p.value")
    print(paste0("logFC > 0: ", nrow(all.results[all.results$p.adjusted < 0.05 & all.results$effect.size > 0,])))
    print(paste0("logFC < 0: ", nrow(all.results[all.results$p.adjusted < 0.05 & all.results$effect.size < 0,])))
    
    print(summary(all.results$effect.size))
  }
}


