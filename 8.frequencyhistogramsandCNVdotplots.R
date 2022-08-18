#This script generates frequency histograms and dot-plots for cn-mops identified CNV locations.
#Last edited on 18/08/22 by Matthew Gaskins.
setwd("D:/correctcnmopsoutput/") #sets working directory containing candidate CNVs window for one chromosome
library("readxl")
tet <- read_excel("tetscaffold14big.xlsx") #manual input of excel file containing cn.mops-generated CNV windows
dip <- read_excel("0.5dipscaff14big.xlsx")
library(dplyr)
library(ggplot2)

tet$middle <- tet$start + (tet$end - tet$start) #identifies middle of each CNV window
dip$middle <- dip$start + (dip$end - dip$start)

tet$ploidy <- 'tetraploid'
dip$ploidy <- 'diploid'
tet$chromno <- 1
dip$chromno <- 1.7 #CHROMOSOME NUMBER NORMALISATION
both <- rbind(tet,dip)

f <- ggplot(both, aes(y=start)) + geom_linerange(aes(xmin=start,xmax=end),linetype=1,color="black") + geom_point(aes(x = start, color=ploidy, size = (end-start))) + geom_point(aes(x=end, color=ploidy, size = (end-start), alpha = mean))
f + scale_alpha(range = c(1, 0.5)) + xlab("CNV Location") + xlim(c(0,30000000))
#plots dotplots of CNVs from each population, with diploids and tetraploids coloured red and turquoise respectively.
#The size and opacity of each CNV point is determined by window width and mean expected log2 fold coverage decrease respectively.
library(ggplot2)


ggplot(both, aes(x=start, fill=ploidy, weights=chromno)) + geom_histogram(alpha=0.5, position="identity") + xlab(c("Start Position (Mb)")) + xlim(c(0,30000000))
#plots frequency histograms comparing CNV locations identified in diploid (red) and tetraploid(turquoise) populations.
