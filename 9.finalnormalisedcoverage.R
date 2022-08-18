#This script generates mean normalized coverage plots for diploid and tetraploid populations from SAMTools-generated read-depth evaluation.
#Last edited by Matthew Gaskins on 18/08/22.
install.packages("readr")
library("readr")
setwd("D:/coverage/coverage2") #sets working directory containing SAMTools-generated file of read-depth at every nucleotide position
dip <- read_table2("alldipsscaff1.coverage") #contains SAMTools-generated read-depth values at every chromosome position of one chromosome
tet <- read_table2("alltetsscaff1.coverage")
library(reshape)
names(dip)[1] <- "scaffold_1" #change name of each column
names(dip)[2] <- "locus"
names(dip)[3] <- "depth"
names(tet)[1] <- "scaffold_1"
names(tet)[2] <- "locus"
names(tet)[3] <- "depth"
#loads the library to rename the column names
library(data.table) # to make sliding window dataframe
library(zoo) # to apply rolling function for sliding window

Xdepth.average<-setDT(dip)[, .(
  window.start = rollapply(locus, width=10000, by=10000, FUN=min, align="left", partial=TRUE),
  window.end = rollapply(locus, width=10000, by=10000, FUN=max, align="left", partial=TRUE),
  coverage = rollapply(depth, width=10000, by=10000, FUN=mean, align="left", partial=TRUE)
), .(scaffold_1)]
#averages read-depth values across 10=kilobase windows to restrict influence of genetic drift
Xdepth.averagetet<-setDT(tet)[, .(
  window.start = rollapply(locus, width=10000, by=10000, FUN=min, align="left", partial=TRUE),
  window.end = rollapply(locus, width=10000, by=10000, FUN=max, align="left", partial=TRUE),
  coverage = rollapply(depth, width=10000, by=10000, FUN=mean, align="left", partial=TRUE)
), .(scaffold_1)]

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
} #calculates mean normalized coverage

Xdepth.averagetet$normalised <- min_max_norm(Xdepth.averagetet$coverage)
Xdepth.averagetet$unlisted <- unlist(Xdepth.averagetet$normalised)
Xdepth.averagetet$middle <- (Xdepth.averagetet$window.start + Xdepth.averagetet$window.end)/2 #calculates middle of each read-depth window

Xdepth.average$normalised <- min_max_norm(Xdepth.average$coverage)
Xdepth.average$unlisted <- unlist(Xdepth.average$normalised)
Xdepth.average$middle <- (Xdepth.average$window.start + Xdepth.average$window.end)/2

plot(Xdepth.average$window.end, Xdepth.average$unlisted, type = 'l', col = "red", xlim = c(25000000,27000000)) #plots diploid mean normalized coverage
lines(Xdepth.averagetet$window.end, Xdepth.averagetet$unlisted, col="blue",lty=2) #adds tetraploid normalized coverage