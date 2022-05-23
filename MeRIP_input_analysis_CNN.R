###################################
#  Cell number-normalized (CNN)   #
#      MeRIP INPUT analysis       #
###################################

#Note: this script is very similar to the CNN RNAseq, but using the ratio human/mouse mapped reads for the normalizing factors (nf)

library(edgeR)

# 1. Formatting tables as needed for downstream analysis

  # Uploading raw tables. This file has been uploaded to the github repository.

Exp <- read.table(file="path")
spike <- c(1,1.325026943,1.230524841,1.882893653,	2.164883885,2.305161538) # human/mouse mapped reads for each sample (relative to sample 1)
N <- colSums(Exp)
d0 <- DGEList(Exp)

  # Filtering low count genes

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
d <- calcNormFactors(d)

  # Adding information about samples for differential analysis

group <- c("FBS", "FBS", "FBS", "Paused", "Paused", "Paused")
batch <- factor(c(1,2,3, 1,2,3))



# 2. Differential gene analysis with edgeR

  # Building the model for limma-voom

mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, lib.size = N * spike, plot=T)
fit <- lmFit(y, mm)
head(coef(fit))

  # Differential analysis for inputs (Paused vs FBS)

contr <- makeContrasts(groupPaused - groupFBS, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

  # Exporting the results

#write.table(y$E, file = "export_folder/MeRIP.input.Counts_log2cpmVoom.CNN.txt")
#write.table(top.table, file = "export_folder/MeRIP.input.toptable.CNN.txt")

