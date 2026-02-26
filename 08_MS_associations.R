# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Load libraries -----------------------------------------------------------
library(phyloseq)
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
library(tidyr)
library(ComplexHeatmap)
library(ggpubr)
library(ggsignif)
library(stringr)

## Local functions

wilcoxon.local = function(data, vars.test, var.group) {
  
  ## Inicialize df results
  res.wilcoxon = data.frame(matrix(data = 0, nrow = 0, ncol = 5))
  colnames(res.wilcoxon) = c("Variable", "Group1", "Group2", "W", "p.value")
  data[,var.group] = factor(data[,var.group])
  combination = levels(unique(data[,var.group]))
  
  for (var in vars.test) {
    df2 = data[!is.na(data[,var]),]
    
    # Extract values for each group
    group1 = df2[df2[,var.group] == combination[1], var]
    group2 = df2[df2[,var.group] == combination[2], var]
    
    # Test de Wilcoxon 
    wilcox_test = wilcox.test(group1, group2, exact = FALSE)
    res.w = wilcox_test$statistic
    res.p = wilcox_test$p.value
    
    # Save results
    res.wilcoxon = rbind(res.wilcoxon, c(var, combination[1], combination[2], res.w, res.p))
  }
  
  colnames(res.wilcoxon) = c("Variable", "Group1", "Group2", "W", "p.value")
  res.wilcoxon$p.adjusted = p.adjust(p = res.wilcoxon$p.value, method = "BH")
  
  res.wilcoxon$p.value = as.numeric(res.wilcoxon$p.value)
  res.wilcoxon$p.adjusted = as.numeric(res.wilcoxon$p.adjusted)
  
  res.wilcoxon[,c(5,6)] = signif(res.wilcoxon[,c(5,6)], digits = 4)
  
  print(DT::datatable(res.wilcoxon, rownames = FALSE, caption = paste0("Wilcoxon test: ", var.group),
                      options = list(pageLength = 30)) %>% 
          formatStyle('p.value', backgroundColor = styleInterval(0.05, c("#c7f9cc", ""))) %>% 
          formatStyle('p.adjusted', backgroundColor = styleInterval(0.05, c("#88d4ab", ""))))
  
  return(res.wilcoxon)
}


kruskal.local = function(data, vars.test, var.group) {
  
  ### Inicialize results
  res.kruskal = data.frame(matrix(data = 0, nrow = 0, ncol = 3))
  colnames(res.kruskal) = c("Variable", "chi-squared", "p.value")
  df = data
  
  ### Kruskal-Wallis & Wilcoxon test
  for (var in vars.test) {
    
    # Test de Kruskal-Wallis para comparar los 3 grupos
    kruskal_test = kruskal.test(df[,var] ~ df[,var.group])
    
    # Save results
    res.c = kruskal_test$statistic
    res.p = kruskal_test$p.value
    
    res.kruskal = rbind(res.kruskal, c(var, res.c, res.p))
  }
  
  
  colnames(res.kruskal) = c("Variable", "chi-squared", "p.value")
  res.kruskal$p.adjusted = p.adjust(p = res.kruskal$p.value, method = "BH")
  
  res.kruskal$`chi-squared` = as.numeric(res.kruskal$`chi-squared`)
  res.kruskal$p.value = as.numeric(res.kruskal$p.value)
  res.kruskal[,c(2:4)] = signif(res.kruskal[,c(2:4)], digits = 4)
  
  print(DT::datatable(res.kruskal, rownames = FALSE, caption = paste0("Kruskal-Wallis test: ", var.group),
                      options = list(pageLength = 30)) %>% 
          formatStyle('p.value', backgroundColor = styleInterval(0.05, c("#c7f9cc", ""))) %>% 
          formatStyle('p.adjusted', backgroundColor = styleInterval(0.05, c("#88d4ab", ""))))
  
  return(res.kruskal)
}


dunn.local = function(data, kruskal, var.group, option, cutoff) {
  
  library(dunn.test)
  
  ## Inicialize results
  res.dunn = data.frame(matrix(data = 0, nrow = 0, ncol = 5))
  colnames(res.dunn) = c("Variable", "Comparisons", "Z", "p.value", "p.adjusted")
  
  if (option == "p.value")  {
    vars = kruskal[kruskal$p.value < cutoff, "Variable"]
  }
  
  if (option == "p.adjusted")  {
    vars = kruskal[kruskal$p.adjusted < cutoff, "Variable"]
  }
  
  
  df = data
  
  ### Dunn test
  for (var in vars) {
    
    # Test de Kruskal-Wallis para comparar los 3 grupos
    dunn_test = dunn.test(x = df[,var], g = df[,var.group], method = "BH")
    
    for (i in 1:length(dunn_test[["comparisons"]])) {
      
      # Save results
      res.c = dunn_test[["comparisons"]][i]
      res.s = dunn_test[["Z"]][i]
      res.p = dunn_test[["P"]][i]
      res.pp = dunn_test[["P.adjusted"]][i]
      
      res.dunn = rbind(res.dunn, c(var, res.c, res.s, res.p, res.pp))
      
    }
  }
  
  colnames(res.dunn) = c("Variable", "Comparisons", "Z", "p.value", "p.adjusted")
  
  res.dunn[,c(3:5)] = apply(res.dunn[,c(3:5)], 2, as.numeric)
  res.dunn[,c(3:5)] = signif(res.dunn[,c(3:5)], digits = 4)
  
  print(DT::datatable(res.dunn, rownames = FALSE,  caption = paste0("Dunn test: ", var.group),
                      options = list(pageLength = 30)) %>% 
          formatStyle('p.value', backgroundColor = styleInterval(0.05, c("#c7f9cc", ""))) %>% 
          formatStyle('p.adjusted', backgroundColor = styleInterval(0.05, c("#88d4ab", ""))))
  
  return(res.dunn)
}


# --- Load data ----------------------------------------------------------------
input.path = "TBD"
PRJEB32762 = readRDS("TBD")

################################################################################
#   Taxa  ~ categorical variables                                              #
################################################################################

taxa.interest = TBD
cat.variables = TBD

metadata = data.frame(PRJEB32762@sam_data)


counts = t(as.data.frame(PRJEB32762@otu_table))

PRJEB32762.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))
data = t(PRJEB32762.vst)[,taxa.interest]

df = merge(metadata, data, by = "row.names")
rownames(df) = df$Row.names

# Filter subset of interest
df = df %>% filter(Condition == "MS")

### STATISTICAL TESTS
tests = taxa.interest

## 2 groups
results.w = wilcoxon.local(data = df, vars.test = tests, var.group = "Treatment")

## +2 groups
results.k = kruskal.local(data = df, vars.test = tests, var.group = "Subtype")
results.d = dunn.local(data = df, kruskal = results.k, var.group = "Subtype", option = "p.adjusted", cutoff = 0.05)

# Parsing to plot Enterotypes
results.d$Group1 = str_split_fixed(string = results.d$Comparisons, pattern = " - ", n = 2)[,1]
results.d$Group2 = str_split_fixed(string = results.d$Comparisons, pattern = " - ", n = 2)[,2]

signif = filter(results.d, p.adjusted < 0.05)

signif$Position = c(14, 12, 14)

signif$Variable = factor(signif$Variable)
signif$Group1 = factor(signif$Group1)
signif$Group2 = factor(signif$Group2)

colnames(signif) = gsub(pattern = "Variable", replacement = "Taxa", x = colnames(signif))

# Plot
df.plot = tidyr::pivot_longer(data = df, cols = taxa.interest,
                              names_to = "Taxa", values_to = "Abundance")

ggplot() + geom_boxplot(data = df.plot, aes(x = Subtype, y = Abundance, fill = Subtype)) +
  theme_bw() + ylab("VST normalized abundance") +
  facet_wrap(~ Taxa, scales = "free", ncol = 3) +
  scale_fill_manual(values = colorss) +
  theme(legend.position = "none", text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + #
  scale_y_continuous(expand = c(0, 0, 0.1, 0)) +  # Mostrar bien el asterisco
  geom_signif(data = signif,
              aes(y_position = Position, xmin = Group1, xmax = Group2, annotations = "*"),
              manual = TRUE, tip_length = 0.01, vjust=0.2, textsize = 4)


################################################################################
#   Taxa  ~ numerical variables                                                #
################################################################################
taxa.interest = TBD
num.variables = TBD

metadata = data.frame(PRJEB32762@sam_data)

counts = t(as.data.frame(PRJEB32762@otu_table))

PRJEB32762.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))
data = t(PRJEB32762.vst)[,taxa.interest]

df = merge(metadata, data, by = "row.names")
rownames(df) = df$Row.names

# FILTER SUBSET OF INTEREST 
df = df %>% filter(Condition == "MS")

g1 = df[, taxa.interest]
g2 = df[,num.variables]

# Correlation test
cor_test_result = cor.test(df$Eisenbergiella, df$Eggerthella, method = "spearman")

results.cor = spearman.local(group1 = g1, group2 = g2, name1 = "Taxa", name2 = "Variables", option = "none")
results.cor$N = as.numeric(results.cor$N)

all.res = correlation.parsing.local(df = results.cor, Group2 = "Variables")

res.rho = as.data.frame(all.res[[1]])
res.rho$Variable = "rho"
res.rho$Taxa = row.names(res.rho)

res.pval = as.data.frame(all.res[[2]])
res.pval$Variable = "p.value"
res.pval$Taxa = row.names(res.pval)

res.padj = as.data.frame(all.res[[3]])
res.padj$Variable = "p.adjusted"
res.padj$Taxa = row.names(res.padj)

all.equal(colnames(res.rho), colnames(res.pval), colnames(res.padj))

  
(p = ggplot(df2plot, aes(x = Variable, y = Taxa, size = abs(rho), colour = Significance)) + 
    geom_point() +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
          text = element_text(size = 16),
          legend.margin = margin(0, 0, 0, 0),
          legend.justification.right = "top",
          plot.margin = margin(t = 40),
          legend.position = "right",
          legend.justification = "center"
    ))

