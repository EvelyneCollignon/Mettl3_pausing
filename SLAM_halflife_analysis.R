############################
#    Half-life analysis    #
#      from SLAMseq        #
############################

#Notes: -Normalized counts were generated with edgeR in the same manner as CNN RNA seq, 
#        but using "T-to-C conversion" raw counts normalized to total raw counts for all genes.


# 1. Formatting tables as needed for downstream analysis

  # Uploading the normalized counts table from SLAMseq. This file has been uploaded to the github repository.
  # Separating WT and Mettl-KO samples in 2 tables and normalizing to time zero


logCPM <- read.table("path", sep='\t')
WT <- log2(2^logCPM[,c(1,4,5,8,11,12)]/2^logCPM$WT_0.1)
KO <- log2(2^logCPM[,c(2,3,6,7,9,10,13,14)]/rowMeans(2^logCPM[,c(2,3)]))

time_WT <- c(0,1,1,3,8,8)
time_KO <- c(0,0,1,1,3,3,8,8)


# 2. Loop to calculate the half-lives separately for WT and KO

  # For WT table

WT_half <- vector()
WT_cor <- vector()

for (i in 1:nrow(WT)) {
  Values <- as.numeric(WT[i,])
  model <- lm(Values ~ 0 + time_WT)
  WT_half[i] <- -1/coef(model)
  WT_cor[i] <- cor(time_WT, as.numeric(WT[i,]))
  if (i%%100 == 0) {print(i)}
}

  # For KO table

KO_half <- vector()
KO_cor <- vector()

for (i in 1:nrow(KO)) {
  Values <- as.numeric(KO[i,])
  model <- lm(Values ~ 0 + time_KO)
  KO_half[i] <- -1/coef(model)
  KO_cor[i] <- cor(time_KO, as.numeric(KO[i,]))
  if (i%%100 == 0) {print(i)}
}

  #Merging WT and KO data

df <- data.frame(
  Gene = rownames(logCPM),
  WT_half = WT_half,
  WT_cor = WT_cor,
  KO_half = KO_half,
  KO_cor = KO_cor
)

  # Removing data that don't make sense or lower confidence in half-life value:
  # -remove negative values (and double zeros) 
  # -remove genes with corr > -0.9

df[df$WT_half < 0.1,'WT_half'] <- 0.1
df[df$KO_half < 0.1,'KO_half'] <- 0.1
df <- df[(df$WT_half < 1 & df$KO_half < 1) == F,]
df <- df[df$WT_cor < -0.8 & df$KO_cor < -0.8,]
df$delta <- df$KO_half - df$WT_half


  # Export results

write.table(df, file="path2", quote=F, sep='\t', row.names = F)



