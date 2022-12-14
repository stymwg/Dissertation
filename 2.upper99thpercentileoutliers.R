##By Matthew Gaskins. Last edited 17/08/22. Extracts top 1% of outliers from ScanTools-generated list of windows.
setwd("D:/") #sets working directory containing list of windows 
all <- read_excel("allscaffolds.xlsx") #Manual Input - requires list of windowed PGM values generated by ScanTools
new_df <- subset(all, Fst_WC<1 & Fst_WC>0) #subsets all FST values between 0 and 1
new_df2 <- subset(new_df, snp_count>25) #extracts all windows with more than 25 SNPs
hist(new_df2$Fst_WC) #plots histogram
d <- density(new_df2$Fst_WC)
plot(d) #plots frequency distribution


cutoff <- round(nrow(new_df2)*0.01) #cutoff for 1% outliers - can be changed depending on desired outlier number
library(dplyr)
arrangedFST <- arrange(new_df2,desc(Fst_WC)) #arranges FST_WC from largest to smallest
FSToutliers_all <- arrangedFST[1:cutoff,]
ggplot(FSToutliers_all, aes(x = window_start, y = Fst_WC, fill = scaffold)) + geom_point(size=2, shape=23) #plots location of identified outliers
FSToutliers <- arrangedFST[1:cutoff,] #extracts 1% outliers
