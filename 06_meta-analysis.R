# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Load libraries -----------------------------------------------------------
library(Biobase)
library(metafor)
library(stringr)
library(limma)


# --- Load data ----------------------------------------------------------------

## General variables
filters = "TBD"
comparison = "TBD"
path = "TBD"

rds.results = list.files(path = path, pattern = paste0("_", comparison, ".rds$"), full.names = TRUE, recursive = TRUE)
individual.results = setNames(lapply(rds.results, readRDS), nm = sapply(rds.results, function(x) strsplit(basename(x), "_")[[1]][1]))


################################################################################
#   Processing                                                                 #
################################################################################
# Taxa
all.taxa = individual.results[[1]]$Taxa
all.taxa = sort(all.taxa, decreasing = FALSE)

# ----------------------- #
# --- logFC matrix ------ #
# ----------------------- #

mat.logFC = matrix(NA, nrow = length(all.taxa), ncol = length(individual.results))
rownames(mat.logFC) = all.taxa
colnames(mat.logFC) = names(individual.results)
head(mat.logFC)

for (i in 1:length(individual.results)) {
  
  dataset = names(individual.results[i])
  df = individual.results[[i]]
  rownames(df) = df$Taxa
  
  ## Order dataframe
  df.ordered = df[all.taxa,]
  print(all.equal(all.taxa, df.ordered$Taxa))
  
  ## Add logFC
  mat.logFC[,dataset] = df.ordered$effect.size
  
}


# ----------------------- #
# --- SD matrix --------- #
# ----------------------- #

mat.SD = matrix(NA, nrow = length(all.taxa), ncol = length(individual.results))
rownames(mat.SD) = all.taxa
colnames(mat.SD) = names(individual.results)
head(mat.SD)

for (i in 1:length(individual.results)) {
  
  dataset = names(individual.results[i])
  df = individual.results[[i]]
  rownames(df) = df$Taxa
  
  ## Order dataframe
  df.ordered = df[all.taxa,]
  print(all.equal(all.taxa, df.ordered$Taxa))
  
  ## Add SD
  mat.SD[,dataset] = abs(df.ordered$SE)
  
}

################################################################################
#   Meta-analysis                                                              #
################################################################################
MA = lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SD[x, ],
                             method = "DL")})

names (MA) = rownames(mat.logFC)
class (MA)
length(MA)
head (MA)
MA[[1]]

################################################################################
#   Explore results                                                            #
################################################################################
cutoff = 0.05

results.MA = as.data.frame(do.call("rbind",
                                    lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub, 
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                            })))

colnames(results.MA) <- c("lower_bound", "logFC", "upper_bound",
                           "p.value", "QE", "QEp", "SE", "tau2", "I2", "H2")

## Corrrect p.value
results.MA$p.adjusted = stats::p.adjust(p = results.MA$p.value, method = "BH")

## Add number of studies where the taxa has been evaluated
sum.na = apply(mat.logFC, 1, is.na)
sum.na = colSums(sum.na==TRUE)

n.studies = ncol(mat.logFC) - sum.na
all.equal(rownames(results.MA), names(n.studies))
results.MA$studies = n.studies
results.MA$Taxa = row.names(results.MA)

results.MA = results.MA[,c(13, 1:12)]

## Significative results
sig.taxa = results.MA[results.MA$p.adjusted < cutoff,]


################################################################################
#   Influence and sensitivity stats                                            #
################################################################################
i = 1

# Add results to the sig.taxa object
for (i in rownames(sig.taxa)){
  
  print(i)
  datasets = colnames(mat.logFC)
  
  # --- INFLUENCE STATS ------------------------------------------------------ #
  sig.taxa[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi) == rep(sign(MA[[i]]$b), length(datasets)))
  
  # how many studies could be influencers?
  inf = influence(MA[[i]])
  res = paste(datasets[inf$is.infl], collapse = ", ")  
  sig.taxa[i, "infl.nstudies"] = ifelse(res =="", "non", res)
  
  # --- SENSITIVITY ANALYSIS ------------------------------------------------- #
  l1 = as.data.frame(leave1out(MA[[i]]))
  rownames(l1) = datasets

  sig.taxa[i, "sensi.global"] <-t.test(x= l1$estimate,
                                           mu=as.numeric(MA[[i]]$b))$p.value
  # number of  studies where pvalue > 0.05 
  res2 <- paste(datasets[l1$pval > cutoff], collapse = ",")  
  sig.taxa[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}


################################################################################
#   Influence and sensitivity plots                                            #
################################################################################

path = "TBD"

while (!is.null(dev.list())) dev.off()

for (mytaxa in rownames(sig.taxa)){
  
  res = rma(yi= mat.logFC[mytaxa,], sei =mat.SD[mytaxa,], method = "DL")
  
  #FOREST PLOT
  forest(res,
         slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
         xlab="Effect size", cex=0.7,
         mlab=paste("DL model for all Studies", sep = " "),
         border = "black", #Color del borde del rombo
         col = "red", #Color del rombo
         main = paste("\n", mytaxa, sep=""))
  text(9,-3, "Effect size [IC 95%]", pos=2, cex = 0.7)

  
  #FUNNEL PLOT
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("Effect size (", mytaxa, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("Effect size (", mytaxa, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("Effect size (", mytaxa, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance",
         back ="darkslategray1",  xlab = paste("Effect size (", mytaxa, ")",sep =""))
  dev.off()
  
  #INFLUENCE PLOTS
  inf = influence(res)
  plot(inf)
}


## Save results

### All results
write.table(x = results.MA, file = paste0(path, "/", comparison, "/all_results_MA.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(object = results.MA, file = paste0(path, "/", comparison, "/all_results_MA.rds"))

### Significant results
write.table(x = sig.taxa, file = paste0(path, "/", comparison, "/significant_results_MA.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(object = sig.taxa, file = paste0(path, "/", comparison, "/significant_results_MA.rds"))
