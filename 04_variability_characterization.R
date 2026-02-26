# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Load libraries --------------------------------------------------------- #
library(ggplot2)
library(phyloseq)
library(vegan)
library(DT)
library(stringr)
library(dplyr)
library(tidyr)
library(dplyr)
library(ggsignif)


## Local functions
pcoa.local= function(ps, method.selected = "PCoA", distance.selected = "bray") {
  
  coord.l = ordinate(ps, method = method.selected, distance = distance.selected)
  
  return(coord.l)
}


capscale.local = function(METADATA, GENUS) {
  
  ## Initialize results
  all = c()
  
  ## Loop single vars dbRDA
  for (i in 1:ncol(METADATA)) { 
    
    ## Analysis
    capsc = vegan::capscale(GENUS ~ METADATA[,i], distance = "bray", na.action=na.exclude)
    an = vegan::anova.cca(capsc, permutations = 9999) 
    
    ## Results
    pval = an["Pr(>F)"][[1]][[1]] 
    Fa = an["F"][[1]][[1]] 
    r2 = vegan::RsquareAdj(capsc)[[1]] 
    r2adj = vegan::RsquareAdj(capsc)[[2]] 
    N = nrow(na.exclude(METADATA[, i, drop=FALSE]))
    all = rbind(all, cbind(Fa, r2, r2adj, N, pval))
  }

  colnames(all) = c("F","r2","r2adj", "N", "p.value") #generate table
  row.names(all) = colnames(METADATA)
  all = as.data.frame(all)

  all$p.adjusted = p.adjust(all$p.value, method="BH")
  
  return(all)
}


capscale.results.local = function(all, names.df) {
  
  all[,c(1:3, 5:6)] = signif(x = all[c(1:3, 5:6)], digits = 3)
  all$ID = rownames(all)
  
  col.order = c("ID", "F", "r2", "r2adj", "N", "p.value", "p.adjusted")
  
  ## Generate table
  DT::datatable(all[,col.order], rownames = FALSE,
                options = list(pageLength = 30)) %>% 
    formatStyle('p.value', backgroundColor = styleInterval(0.05, c("#c7f9cc", ""))) %>% 
    formatStyle('p.adjusted', backgroundColor = styleInterval(0.05, c("#88d4ab", "")))
}


capscale.stepwise.local = function(all, METADATA, GENUS, PTRESH = 0.05, diss_metric = "bray") {

  sign_cap = row.names(all[which(all$p.value < PTRESH),])
  paste0(colnames((METADATA[,sign_cap])), collapse = ", ")
  
  MET = data.frame(METADATA[,sign_cap]) 
  
  MET = na.exclude(MET) 
  GEN = GENUS[which(row.names(GENUS) %in% row.names(MET)), ]
  
  distmat = vegdist(GEN, method=diss_metric)
  
  ## Capsale models
  mod0 = capscale(distmat ~ 1) #H0: unconstrained ordination
  mod1 = capscale(distmat ~ ., data = MET) #H1: full constrained ordination, all metadata in MET
  
  attach(MET)
  
  ## Stepwise results
  step.res = ordiR2step(mod0, scope=formula(mod1), data=MET, direction="forward", Pin = 0.1, R2scope = TRUE, pstep = 100, perm.max = 999, permutations=999, trace = T) #forward stepwise dbRDA
  
  ## Parse results
  res = step.res$anova
  row.names(res) = gsub(pattern="\\+ ", "", row.names(res))
  colnames(res) = paste0("RDAcumul_",colnames(res))
  
  return(res)
}


capscale.stepwise.results.local = function(all, res, MET, names.df) {
  
  all$ID = rownames(all)

  all = data.frame(merge(all,res,by="row.names",all=T),row.names=1)
  
  ## Parse complete results
  all = all[order(all$r2,decreasing=TRUE),]
  all = all[order(all$RDAcumul_R2.adj),]
  all[,"RDAcumul_N"] = nrow(MET)
  
  toround = c("F", "r2", "r2adj", "p.adjusted", "RDAcumul_R2.adj", "RDAcumul_AIC", "RDAcumul_F")
  all[,toround] = signif(x = all[,toround], digits = 3)
  
  all[is.na(all$ID), "ID"] = "<All variables>"
  
  col.order = c("ID", "F", "r2", "r2adj", "N", "p.value", "p.adjusted", "RDAcumul_R2.adj", "RDAcumul_Df", "RDAcumul_AIC",
                "RDAcumul_F", "RDAcumul_Pr..F.", "RDAcumul_N")
  
    print(DT::datatable(all[,col.order], rownames = FALSE,
                      options = list(pageLength = 30)) %>% 
          formatStyle('p.value', backgroundColor = styleInterval(0.05, c("#c7f9cc", ""))) %>% 
          formatStyle('p.adjusted', backgroundColor = styleInterval(0.05, c("#88d4ab", ""))) %>% 
          formatStyle('RDAcumul_Pr..F.', backgroundColor = styleInterval(0.05, c("#c7f9cc", "")))) 
  
  return(all)
}

# --- Load data ----------------------------------------------------------------
ps = readRDS("TBD")

################################################################################
#   Alpha diversity                                                            #
################################################################################

# Shannon diversity
colorss = c("TBD")
alpha.tab = estimate_richness(ps, measures = "Shannon")
df = merge(ps@sam_data, alpha.tab, by = "row.names")

# Test 
results.k = kruskal.local(data = df, vars.test = "Shannon", var.group = "TBD")
results.d = dunn.local(data = df, kruskal = results.k, var.group = "TBD", option = "p.value", cutoff = 0.25)

# Plot
ggplot() + geom_boxplot(data = df, mapping = aes(x = TBD, y = Shannon, fill = TBD)) +
  theme_bw() +
  scale_fill_manual(values = colorss) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  geom_signif(data = signif,
              aes(y_position = Position, xmin = Group1, xmax = Group2, annotations = "*"),
              manual = TRUE, tip_length = 0.01, vjust=0.2, fontface = "bold")


################################################################################
#   Beta diversity                                                             #
################################################################################
data = ps
coords = pcoa.local(ps = data)

# Plot
plot_ordination(data, coords, color = "TBD") +
    geom_point(size = 1) +
    stat_ellipse(type = "t", linetype = 3, geom = "polygon", aes(fill = "TBD"),  alpha = 0.01) +
    theme_bw() +
    tune::coord_obs_pred() + # Same X and Y scale
    scale_color_manual(values = colorss) +
    scale_fill_manual(values = colorss) +
    theme(text = element_text(size = 16))



################################################################################
#   dbRDA                                                                      #
################################################################################
colorss = c("TBD")

# --- Univariate ---------------------------------------------------------------

# Metadata selection
metadata2explore = c("TBD")
metadata = data.frame(data@sam_data)
META = metadata[,metadata2explore]

## Change unknown to NA
for (i in 1:ncol(META)) {
  x = which(META[,i] == "unknown")
  META[x, i] = NA
}

# Genera count matrix
GENERA = otu_table(data) 

# Check correct order of names
GENERA = GENERA[rownames(META),]
META = META[rownames(GENERA),]
all.equal(rownames(META), rownames(GENERA))

# dbRDA
results.capscale = capscale.local(METADATA = META, GENUS = GENERA)
capscale.results.local(all = results.capscale)


## Plot
results.capscale.plot = results.capscale
results.capscale.plot$Variable = rownames(results.capscale.plot)

### Order & factor
results.capscale.plot = results.capscale.plot[order(results.capscale.plot$r2adj, decreasing = TRUE),]
results.capscale.plot$Variable = factor(results.capscale.plot$Variable, levels = rev(results.capscale.plot$Variable))

### Plot
ggplot(data = results.capscale.plot) + geom_col(mapping = aes(x = Variable, y = r2adj*100), fill = colorss[1]) +
  coord_flip() +
  geom_text(aes(x = Variable, y = r2adj*100, 
                label = ifelse(p.value < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.2, size = 7) +
  theme_bw() +
  ylab("Univariate variance (%)") +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 12)) +
  expand_limits(y = max(results.capscale.plot$r2adj * 100) * 1.05)

# --- Cumulative ---------------------------------------------------------------
results.stepwise = capscale.stepwise.local(all = results.capscale, METADATA = META, GENUS = GENERA, PTRESH = 0.05, diss_metric = "bray")
all.results = capscale.stepwise.results.local(all = results.capscale, res = results.stepwise, MET = META, names.df = NULL)


### Plot univariate + cumulative
all.results$Variable = rownames(all.results)

df2plot = bind_rows(
  all.results %>%
    select(R2 = r2adj, p.value = p.value, Variable) %>%
    mutate(Variance = "Univariate"),
  
  all.results %>%
    select(R2 = RDAcumul_R2.adj, p.value = RDAcumul_Pr..F., Variable) %>%
    mutate(Variance = "Cumulative")
) %>%
  select(Variance, R2, p.value, Variable) %>%
  filter(complete.cases(.))

df2plot = as.data.frame(df2plot)

df2plot$Variance = factor(df2plot$Variance, levels = rev(c("Univariate", "Cumulative")))

df2plot$R2 = ifelse(df2plot$R2 < 0, 0, df2plot$R2)

max.cumul = max(df2plot$R2[df2plot$Variance == "Cumulative"], na.rm = TRUE)
vars_with_cumul = unique(df2plot$Variable[df2plot$Variance == "Cumulative"])
vars_without_cumul = setdiff(unique(df2plot$Variable), vars_with_cumul)

## Cumulative shadow
shadows.cumul = data.frame(
  Variance = "Cumulative_shadow",
  R2 = max.cumul,
  p.value = 1,
  Variable = vars_without_cumul
)

df2plot = bind_rows(df2plot, shadows.cumul)

# Order variables
original_cumul_order = df2plot %>%
  filter(Variance == "Cumulative") %>%
  arrange(R2) %>%
  pull(Variable)

vars_final_order = c(original_cumul_order, setdiff(unique(df2plot$Variable), original_cumul_order))
df2plot$Variable = factor(df2plot$Variable, levels = rev(vars_final_order))

## Order
ggplot(df2plot) +
  geom_col(aes(x = R2 * 100, y = Variable, fill = Variance),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  
  # Línea de referencia vertical
  geom_vline(xintercept = max.cumul*100, linetype = "dashed", color = "black") +
  
  # Asteriscos de significancia
  geom_text(aes(x = R2 * 100, y = Variable, group = Variance,
                label = ifelse(!is.na(p.value) & p.value < 0.05 & Variance != "Cumulative_shadow", "*", "")),
            position = position_dodge2(width = 0.9), vjust = 0.5, hjust = -0.2, size = 5) +
  
  # Colores y leyenda
  scale_fill_manual(
    values = c(
      "Univariate" = "#b95f89",
      "Cumulative" = "#01497c",
      "Cumulative_shadow" = "#d7e3fc"
    ),
    breaks = c("Univariate", "Cumulative") # Ocultar shadow de la leyenda
  ) +
  theme_bw() +
  expand_limits(x = max(df2plot$R2 * 100) * 1.05) +
  xlab("Variance (%)") + ylab("Variables") +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12)
  )



