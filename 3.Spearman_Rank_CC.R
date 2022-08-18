##Calculates Spearman's Rank Correlation Coefficient between PGMs and SNP count within windows.
#Last edited by Matthew Gaskins on 17/08/22.
install.packages("readr") #installs readr package
library("readr")
setwd("D:/") #set working directory to that containing list of all ScanTools-generated window
library(readxl)
FST <- read_excel("allscaffolds.xlsx") #reads in list of ScanTools generated windows from Excel
library("dplyr")
plot(FST$snp_count, FST$Fst_WC) #plots distribution of outliers as Fst versus SNP count
corr <- cor.test(x=FST$Fst_WC, y=FST$snp_count, method = 'spearman') #calculates Spearman's Rank Correlation Coefficient
corr
abline(lm(FST$Fst_WC ~ FST$snp_count), col = "red", lwd = 3)










new_df2 <- subset(new_df, N_VARIANTS>25)
hist(new_df2$WEIGHTED_FST)
d <- density(new_df2$WEIGHTED_FST) #change y lab to frequency
plot(d)


cutoff <- round(nrow(new_df2)*0.01)

arrangedFST <- arrange(new_df2,desc(WEIGHTED_FST))
FSToutliers_SOL <- arrangedFST[1:cutoff,]






plot(ctg1$BIN_START, ctg1$WEIGHTED_FST) #for all contigs






install.packages("writexl")
library("writexl")
write_xlsx(FSToutliers_SOL,"FST_PAULLAN_outliers_vcftools.xlsx")
