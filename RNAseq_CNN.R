###################################
#  Cell number-normalized (CNN)   #
#  RNAseq analysis (using ERCCs)  #
###################################

library(edgeR)


# 1. Formatting tables as needed for downstream analysis

  # Uploading raw tables (for genes and ERCC spikes). These files have been uploaded to the github repository.

Exp <- read.table("path", row.names = 1)
ercc_expr <- read.table("path", row.names = 1)
d0 <- DGEList(Exp)

  # Filtering low count genes

cutoff <- 1           
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
d <- calcNormFactors(d)
N <- colSums(Exp)
nf <- calcNormFactors(ercc_expr, lib.size=N) #calculating normalizing factors (nf) based on ERCCs
dim(d) #number of genes left

  # Adding information about samples for differential analysis

ttmt <- c(rep(c("FBS", "FBS", "Paused","Paused"), 3))
batch <- c(rep(1, 4), rep(2, 4), rep(3, 4))
genotype <- c(rep(c("WT", "KO"), 6))
group <- interaction(ttmt, genotype)



# 2. Differential gene analysis with edgeR

  # Building the model for limma-voom

mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, lib.size = N * nf, plot = T)  #lib.size adjusted by nf for CNN purposes
fit <- lmFit(y, mm)
head(coef(fit))

  # Differential analysis for Mettl3 +/+ cells

contr <- makeContrasts(groupPaused.WT - groupFBS.WT, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
WTtop.table <- top.table
WTtop.table$GeneID <- row.names(WTtop.table)

  # Differential analysis for Mettl3 -/- cells
contr <- makeContrasts(groupPaused.KO - groupFBS.KO, levels = colnames(coef(fit))) #Change if other comparison
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many are significant?
KOtop.table <- top.table
KOtop.table$GeneID <- row.names(KOtop.table)

logCPM <- data.frame(y$E) #saving normalized counts for all samples
top.table <- full_join(WTtop.table, KOtop.table, by="GeneID") #merging the 2 top tables

  # Exporting the results

#write.table(logCPM, file = "export_folder/Counts_log2cpmVoom.CNN.txt")
#write.table(top.table, file = "export_folder/toptable.CNN.txt")

