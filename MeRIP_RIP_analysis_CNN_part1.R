###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part1    #
###################################

#Notes: -Raw counts were extracted for all peaks from bam files with Diffbind package.
#       -Part1 normalizes peaks values in a CNN approach.
#       -Part2 (other script) adjusts RIP values to consider changes in genes expression (based on input CNN values)

library(edgeR)
library(tidyverse)

# 1. Formatting tables as needed for downstream analysis

  # Uploading raw tables for peaks. This file has been uploaded to the github repository.
  # The peak information (seqnames, start, end, GeneID) is saved separately and will be merged back later

Peaks <- read.table(file="path")
Ref <- Peaks[,c(1:3,10)]
rownames(Peaks)<- paste("ref", 1:nrow(Ref))
Ref$ID <- paste("ref", 1:nrow(Ref))
Peaks <- unique(Peaks[,-10]) # removing duplicated peaks, i.e. those that mapped to more than 1 gene
spike <- c(1,	1.33564167,	1.237347774,	1.53725346,	1.944077862,	1.998910646) # human/mouse mapped reads for each sample (relative to sample1)
d0 <- DGEList(Peaks[,c(7:9,4:6)]) #reordering columns (FBS first, paused second)
N <- colSums(Peaks[,c(7:9,4:6)])  #reordering columns (FBS first, paused second)

  # Filtering low count peaks

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

  # Bring back peak information and export

counts <- data.frame(y$E)
counts$ID <- rownames(counts)
counts <- inner_join(counts, Ref, by=c("ID"))
#write.table(counts, file = "export_folder/MeRIP.RIP.Counts_log2cpmVoom.CNN.txt")

top.table$ID <- rownames(top.table)
top.table <- inner_join(top.table, Ref, by=c("ID"))
#write.table(top.table, file = "export_folder/MeRIP.RIP.toptable.CNN.txt")

