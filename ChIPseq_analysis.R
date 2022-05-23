###################################
#  Cell number-normalized (CNN)   #
#      Mettl3 ChIP analysis       #
###################################

#Notes: -Peak cpm counts were extracted for from bam files with Diffbind package. This script performs cell number-normalization.
#       -TSS analysis was performed in the same manner but using a 1kb window surround the most upstream TSS for all Refseq-annotated genes.

library(tidyverse)
library(edgeR)

# 1. Formatting tables as needed for downstream analysis

  # Uploading count tables for peaks. This file has been uploaded to the github repository.
  # The peak information (seqnames, start, end, GeneID) will be merged back later. 

counts <- read.table("path" )
d0 <- DGEList(counts[,6:9])
spike <- c(0.579, 0.841, 0.556, 0.739) # human/mouse mapped reads for each sample (relative to respective input sample)

 # Filtering low count peaks

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
d <- calcNormFactors(d)
dim(d) # number of peaks left

 # Adding information about samples for differential analysis

ttmt <- c("FBS","FBS", "Paused", "Paused")
batch <- c(1,2,1,2)
group <- interaction(ttmt)



# 2. Differential gene analysis with edgeR

 # Building the model for limma-voom

mm <- model.matrix(~0 + group)
y <- voom(d, mm, lib.size = 10^6 * spike, plot = F)
fit <- lmFit(y, mm)
head(coef(fit))

 # Differential analysis for inputs (Paused vs FBS)

contr <- makeContrasts(groupPaused - groupFBS, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

 # Bring back peak information and export

top.table$ID <- as.numeric(rownames(top.table))
top.table <- left_join(top.table, counts[,c(1:3,10,12)], by=c("ID"))
norm <- data.frame(y$E) %>% mutate(ID = as.numeric(rownames(y$E))) # adding individual values to the top table 
                                                                   # M3.1 and 2 are FBS, M3.3 and 4 are paused
top.table <- top.table %>% left_join(norm, by=c("ID"))
#write.table(top.table, "export_path/ChIP_results.txt")
