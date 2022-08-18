#This script calculates the significance of overlap between candidates of differentiation between multiple populations
#Last edited by Matthew Gaskins on 18/08/22.
BiocManager::install("GeneOverlap")
library("GeneOverlap")
setwd("D:/") #sets working directory to that containing list of candidate TAIR IDs from multiple populations
library(readr)
outliers <- read_excel("externalcandidates.xlsx") #manual input - an excel file containing columns of candidate TAIR IDs from multiple populations
go.obj <- newGeneOverlap(outliers$arenosa,outliers$internal)
go.obj <- testGeneOverlap(go.obj)
go.obj #outputs significance of overlap value
