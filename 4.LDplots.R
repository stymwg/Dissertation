##Plots LD across all scaffolds in different colours.
#Last edited by Matthew Gaskins on 17/08/22.
setwd("D:/") #sets working directory to that containing long .LD file
install.packages("readr")
library(readr)
linkage <- read_table2("out.hap.ld") #contains LD values for all scaffolds in one long column in the format output by VCFTools. Another column ('CHR') provides scaffold name. 
linkage[is.na(linkage)] <- 0 #removes NA values
library(ggplot2) #may require pre-installation
library(dplyr)
install.packages("devtools")
library(devtools)
devtools::install_github("kassambara/easyGgplot2", force = TRUE)
library(easyGgplot2)
table <- as.data.frame(linkage)
table$Distance <- table$POS2 - table$POS1 #locates middle position of window
install.packages("reshape")
library(reshape)
table<-rename(table,c("R^2"="Rsquared")) #changes name of this column to Rsquared

ggplot2.scatterplot(data=table, xName="Distance",yName="Rsquared", 
                    size=3, groupName="CHR", groupColors=c('black','pink','orange', 'yellow', 'turquoise', 'red', 'blue','magenta','brown','grey','cyan','purple','light green','dark green'), addRegLine=TRUE, smoothingMethod = "loess", regLineSize = 1.5) #plots LD for each contig on same graph
#plots LD across each scaffold
