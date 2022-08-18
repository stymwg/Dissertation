##Generates absolute AFD for candidate windows.
#Last edited by Matthew Gaskins on 18/08/22,
setwd("D:/")
library(readr)
rm(list=ls())
# to start, set up a new environment by calling some packages:
library(ggplot2)
library(tidyverse)

# create this table for each population using GATK
pop1 <- read_table2("outputdip.table") #manual input
pop2 <- read_table2("outputtet.table")
library(dplyr)
results1 <- pop1 %>% filter(CHROM == "scaffold_12") #extracts windows from one scaffold
results2 <- pop2 %>% filter(CHROM == "scaffold_12")
DIP <- c("CHROM", "POS", "14-QUA-54-3.AD", "14-QUA-54-3.DP", "15-QUA-3-1.AD", "15-QUA-3-1.DP", "15-QUA-5-3.DP", "15-QUA-5-3.AD", "15-QUA-7-1.AD", "15-QUA-7-1.DP", "15-QUA-35-3.DP", "15-QUA-35-3.AD","15-QUA-42-1.AD","15-QUA-42-1.DP","15-QUA-51-5.AD","15-QUA-51-5.DP")
TET <- c("CHROM", "POS", "15-QUA-1-3.AD", "15-QUA-1-3.DP", "15-QUA-28-1.AD", "15-QUA-28-1.DP", "15-QUA-30-4.DP", "15-QUA-30-4.AD", "15-QUA-22-2.AD", "15-QUA-22-2.DP", "15-QUA-26-4.DP", "15-QUA-26-4.AD","15-QUA-49-5.AD","15-QUA-49-5.DP")
DIP_table <- results1[,DIP] #only retains column relevant to each population
TET_table <- results2[,TET]

DIP_table$'corrected_14-QUA-54-3.AD' <- as.numeric(sub(".*,", "", DIP_table$'14-QUA-54-3.AD')) #each line of this region of code creates a new column containing the values after the comma from the AD column for each individual DIP
DIP_table$'corrected_15-QUA-3-1.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-3-1.AD'))
DIP_table$'corrected_15-QUA-5-3.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-5-3.AD'))
DIP_table$'corrected_15-QUA-7-1.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-7-1.AD'))
DIP_table$'corrected_15-QUA-35-3.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-35-3.AD'))
DIP_table$'corrected_15-QUA-42-1.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-42-1.AD'))
DIP_table$'corrected_15-QUA-51-5.AD' <- as.numeric(sub(".*,", "", DIP_table$'15-QUA-51-5.AD'))



TET_table$'corrected_15-QUA-1-3.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-1-3.AD')) #each line of this region of code creates a new column containing the values after the comma from the AD column for each individual DIP
TET_table$'corrected_15-QUA-28-1.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-28-1.AD'))
TET_table$'corrected_15-QUA-30-4.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-30-4.AD'))
TET_table$'corrected_15-QUA-22-2.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-22-2.AD'))
TET_table$'corrected_15-QUA-26-4.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-26-4.AD'))
TET_table$'corrected_15-QUA-49-5.AD' <- as.numeric(sub(".*,", "", TET_table$'15-QUA-49-5.AD'))


DIP_table$total_AD <- DIP_table$'corrected_14-QUA-54-3.AD' + DIP_table$'corrected_15-QUA-3-1.AD' + DIP_table$'corrected_15-QUA-5-3.AD' + DIP_table$'corrected_15-QUA-7-1.AD' + DIP_table$'corrected_15-QUA-35-3.AD' + DIP_table$`corrected_15-QUA-42-1.AD` + DIP_table$`corrected_15-QUA-51-5.AD`
DIP_table$total_DP <- DIP_table$'14-QUA-54-3.DP' + DIP_table$'15-QUA-3-1.DP' + DIP_table$'15-QUA-5-3.DP' + DIP_table$'15-QUA-7-1.DP' + DIP_table$'15-QUA-35-3.DP' + DIP_table$`15-QUA-42-1.DP` + DIP_table$'15-QUA-51-5.DP'
DIP_table$DIPAF <- DIP_table$total_AD / DIP_table$total_DP #creates total DIP allele frequency values by dividing total AD values from all DIP populations by the total DP values from all DIP populations

TET_table$total_AD <- TET_table$'corrected_15-QUA-1-3.AD' + TET_table$'corrected_15-QUA-28-1.AD' + TET_table$'corrected_15-QUA-30-4.AD' + TET_table$'corrected_15-QUA-22-2.AD' + TET_table$'corrected_15-QUA-26-4.AD' + TET_table$'corrected_15-QUA-49-5.AD'
TET_table$total_DP <- TET_table$'15-QUA-1-3.DP' + TET_table$'15-QUA-28-1.DP' + TET_table$'15-QUA-30-4.DP' + TET_table$'15-QUA-22-2.DP' + TET_table$'15-QUA-26-4.DP' + TET_table$`15-QUA-49-5.DP`
TET_table$TETAF <- TET_table$total_AD / TET_table$total_DP #creates total TET allele frequency values by dividing total AD values from all TET populations by the total DP values from all TET populations
TET_table$AFD <- abs(TET_table$TETAF - DIP_table$DIPAF) #outputs the positive difference between allele frequencies of diploid and tetraploid populations at each scaffold position

plot(TET_table$POS, TET_table$AFD, ylab = "Absolute AFD", xlab = "Chromosome Position (Mb)", xlim=c(17745000,17750000), xaxs="i") #plots the absolute allele frequency differences between diploid and tetraploid populations

