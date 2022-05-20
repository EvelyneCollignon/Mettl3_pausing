###################################
#  Cell number-normalized (CNN)   #
#   MeRIP RIP analysis - Part2    #
###################################

#Note: The results from this script are considered as the fully normalized m6A values, adjusted for gene expression changes and based on CNN approach.

library(tidyverse)
library(edgeR)


# 1. Formatting tables as needed for downstream analysis

  # Uploading cell number-normalized counts from input and RIP peak analyses.

  # RIP and Input files results from the previous MeRIP scripts (MeRIP RIP Part1 and MeRIP input scripts). 
  # Annotation of the peaks was manually verified after running the "Part1" script.
  # Starting files are provided in the repository for reproducibility (Input_CNN.txt and RIP_CNN.txt).


Input <- read.table("path")
Input$GeneID <- rownames(Input)
RIP <- read.table("path")
Norm.RIP <- inner_join(RIP, Input, by="GeneID")
Norm.RIP$Ave.input <- log2(rowMeans(2^Norm.RIP[,12:17]))

  # Adjusting RIP values for gene expression
  # All RIP values are corrected by dividing by the ratio (Input.sample/Input.average) (note: all values are in log2 format so this appears as "-Input#+Ave.Input")

Norm.RIP <- mutate(Norm.RIP, 
                   Ratio1 = RIP1-Input1+Ave.input,
                   Ratio2 = RIP2-Input2+Ave.input,
                   Ratio3 = RIP3-Input3+Ave.input,
                   Ratio4 = RIP4-Input4+Ave.input,
                   Ratio5 = RIP5-Input5+Ave.input,
                   Ratio6 = RIP6-Input6+Ave.input)
  
  # Formatting the "Ratio" table (i.e. the adjusted RIP peaks used for the final analysis)
  # No filtering is performed as both input and RIP data are already filtered for low counts

Ratio <- cbind(Norm.RIP[,c("seqnames", "start", "end", "GeneID")], 2^(Norm.RIP[,19:24]) )
Ratio$RefID <- paste("Ref", 1:nrow(Ratio)) # used to bring back peak information later
rownames(Ratio) <- Ratio$RefID
d0 <- DGEList(Ratio[,5:10])
d <- calcNormFactors(d0)
N <- colSums(Ratio[,5:10])

  # Adding information about samples for differential analysis

group <- c("FBS", "FBS", "FBS", "Paused", "Paused", "Paused")
batch <- factor(c(1,2,3, 1,2,3))



# 2. Differential gene analysis with edgeR

  # Building the model for limma-voom
  # As the data is already fully normalized, further re-scaling is prevented by fixing the library size in the voom function

mm <- model.matrix(~0 + group + batch)
y <- voom(d, mm, plot = T, lib.size = rep(10^6,6))  #lib.size fixed as 10^6 for all samples
fit <- lmFit(y, mm)
head(coef(fit))

  # Differential analysis for inputs (Paused vs FBS)

contr <- makeContrasts(groupPaused - groupFBS, levels = colnames(coef(fit))) 
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > log2(1.5))) #how many peaks are significant?

  # Bring back peak information and export

top.table$RefID <- rownames(top.table)
top.table <- inner_join(top.table, Ratio, by=c("RefID"))
top.table <- top.table[,-c(7,18)]
#write.table(top.table, file = "export_folder/m6A.CNN.toptable.txt")

All <- inner_join(Norm.RIP, top.table[,1:10], by=c("GeneID", "seqnames", "start", "end"))
All <- All[,c(8:11, 25:29, 12:17, 1:6, 19:24)]
#write.table(All, file = "export_folder/m6A.CNN.all_results.txt")
