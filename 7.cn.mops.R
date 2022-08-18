### R code from vignette source 'cn.mops.Rnw'
### Generates list of candidate CNVs from one scaffold at a time
#Last edited by Matthew Gaskins on 18/08/22.
###################################################
### code chunk number 1: cn.mops.Rnw:41-48
###################################################
options(width=75)
set.seed(0)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("cn.mops") #installs cn.mops package
library(cn.mops)
library(Biobase)
library(GenomicRanges)
library(GenomeInfoDb)
cn.mopsVersion <- packageDescription("cn.mops")$Version

setwd("D:/alldips1tet") #working directory containing all diploid bams. For tetraploid analysis, one tetraploid is also present in this directory.
BAMFiles <- list.files(pattern=".bam$") #Creates list of bam files in the directory. Needs multiple BAM files
bamDataRanges <- getReadCountsFromBAM(BAMFiles, refSeqNames="scaffold_7", WL=2000) #identifies normalized read count changes in 2000+ base windows
res <- cn.mops(bamDataRanges,lowerThreshold = -0.1) #change lower threshold to reflect ploidy
res
resCNMOPS <- calcIntegerCopyNumbers(res) #calculates copy number changes in each identified coverage drop
resCNMOPS
cnvs(resCNMOPS) #list of CNVs
cnvr(resCNMOPS) #list of common variant regions
ranges(cnvr(resCNMOPS))
segplot(resCNMOPS) #plots position of detected CNVs across one chromosome for each individual
plot(resCNMOPS,which=2) #plots normalized read count ratios between individuals
print(resCNMOPS)

results.CNVs  <- as.data.frame(cnvs(resCNMOPS)) #converts identified copy number changes to a data frame
library("writexl")
write_xlsx(results.CNVs,"scaffold7_WL500.xlsx") #writes candidate CNVs to an excel file


