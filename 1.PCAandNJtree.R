#===============================================================================
#                               Welcome!
# This adegenet script was provided by Levi Yant and altered by our lab group.
# From an initial VCF file (from GATK) it does a PCA, K-means clustering
# and also calculates a matrix of genetic distances and does an AMOVA.
#
# It was written for plants so it can deal with different ploidies (useful!)
#
# If a comment doesn't have name, it was made by Levi, otherwise it has the
# name of the person that has done it between hashtags
#
# Script by: Filip Kolar 2017, further edits by Sian Bray 2018 & Levi Yant 2022
# altered by: Ana C. da Silva and Matthew Gaskins
# Date: May 2022
#===============================================================================

#Ana# set working directory if different than what you are using
#setwd("/Users/leviyant/Dropbox/Teaching/Semester2/project/DATA/test")

#install packages if you don't have them!
install.packages("adegenet", dep=TRUE)
install.packages("StAMPP")

#Ana# this setting should print warnings as they occur
options(warn=1)

#Ana# call the libraries needed:
library(adegenet)
library(StAMPP)
install.packages("vcfR")
library(vcfR)
library(ggplot2)
library(MASS)
library(adegraphics) #not strictly necessary for all of this (homebrew R installs will interfere)
#library(pegas) #Ana# Loaded automatically with StAMPP
#library(ape) #Ana# Loaded automatically with StAMPP
#library(ade4) #Ana# Loaded automatically with adegenet!

######################=========MODIFIED FUNCTIONS=========######################

# a function for conversion from vcfR object to genlight in tetraploids
##Levi##: note not all of this is necessary for LIFE4136 project, but some is helpful
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0     #Ana# for diploids it's only lines with hashtag below
  x[x == "0|1"] <- 1     #diploid#
  x[x == "1|0"] <- 1     #diploid#
  x[x == "1|1"] <- 2     #diploid#
  x[x == "0/0"] <- 0     #diploid#
  x[x == "0/1"] <- 1     #diploid#
  x[x == "1/0"] <- 1     #diploid#
  x[x == "1/1"] <- 2     #diploid#
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}
# a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks...
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS to support thousands of samples,
  # this could be replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  
  ## PERFORM THE ANALYSIS ## ---------------------------------------------------
  # eigen analysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  # scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  # rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  
  ## GET LOADINGS ## -----------------------------------------------------------
  # need to decompose X^TDV into a sum of n matrices of dim p*r
  # but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ## ----------------------------------------------------------
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
######################====================================######################

# IMPORT SNP data from VCF
#vcf <- read.vcfR("BF_65mb_F4.4dg_no_new_PA_5.vcf.gz", nrows=10000)      # nrow = number of SNPs to read in
# vcf <- read.vcfR("BF_65mb_F4.4dg_no_new_PA_5.vcf.gz")   #read in all data

# here Levi used the .gz file, but we are going to use the vcf directly
setwd("D:/vcfs/")
vcf <- read.vcfR("finalcorrectrenamedreheadered.vcf")  #Ana# Manual Input!!!
# if you want to be crazy, you can import the full dataset!
#vcf <- read.vcfR("Ana_renamedfull_noPA5.vcf")

# convert to genlight 	
##Levi## this uses the modified function vcfR2genlight.tetra (see Modified functions section)
aa.genlight <- vcfR2genlight.tetra(vcf)
#Ana# here the bit "tetra" just means that this could be used with tetraploids

locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")  # add real SNP.names
## ----
# Ana explains:
# We can use the feature below to group the populations by their characteristics!
# However to do that, I had to rename the populations (on the vcf) so that the
# first two characters now mean inland (IN) or coastal (CO), and then there are 
# three characters that identify location, and one character that identifies sample nr
# So if you choose 2, you get Inland vs Coastal, if 5 grouping by location,
#and 6 for all individuals (no grouping) - which is very useful!
## ----
pop(aa.genlight) <-substr(indNames(aa.genlight),1,3)  # add pop names: here pop names are first chars of ind name

#check    =====VERY IMPORTANT===
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)

######################====================================######################
#   PCA     --------------------------------------------------------------------

#Matt# added this to get rid of the missing values:
toRemove <- is.na(glMean(aa.genlight, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove) # position of entirely non-typed loci
aa.genlight_correct <- aa.genlight[, !toRemove]

#run PCA
# this uses the modified function glPcaFast to run faster!
pca.1 <- glPcaFast(aa.genlight_correct, nf=300)

#Plotting PCA
scatter(pca.1, posi="topleft")  # plot scatter with the individual labels
title("PCA of the population's comparison")
loadingplot(pca.1)

# proportion of explained variance by each axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis
pca.1$eig[4]/sum(pca.1$eig) # proportion of variation explained by 4th axis
pca.1$eig[5]/sum(pca.1$eig) # proportion of variation explained by 5th axis
pca.1$eig[6]/sum(pca.1$eig) # proportion of variation explained by 6th axis
pca.1$eig[7]/sum(pca.1$eig) # proportion of variation explained by 7th axis

#Ana# you can adjust "xax" and "yax" for the PCAs you want to compare (here use just 1 and 2)
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=myCol2)

myCol2=alpha(c('green', 'red', 'blue')) #sets colours of populations


# save nice figs
pdf ("PCA_populations.pdf", width=14, height=7)
#Ana# I've been trying to move that data label from the middle of the group /(T_T)\ HELP!
g1 <- s.class(pca.1$scores, pop(aa.genlight), xax=1, yax=2, col=myCol2,
                plabels.box.draw =F, plabels.cex = 0.9, paxes.draw=T,
                ylab="PC2", xlab="PC1", ppoints.cex=0.5, ppoints.col=myCol2)
#g1 #Ana# I like to see the graphs asap - but comment out before doing the pdf
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = myCol2,
               plabels = list(box = list(draw = FALSE), optim = TRUE),
               ylab="PCA1", xlab="PCA2", paxes.draw=T, pgrid.draw=F, plabels.cex=0.8, plot = FALSE)
#g2 #Ana# I like to see the graphs asap - but comment out before doing the pdf
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#===============================================================================

#Ana# I've used the tutorial "Analysing genome-wide SNP data using adegenet 2.1.5"
# to do more detailed graphs, as described in pages 42/43/44/45
#===============================================================================
#first we start by doing a NJ tree
NJtree <- nj(dist(as.matrix(aa.genlight_correct)))
NJtree
# I've pasted here the output:
## 
## Unrooted; includes branch lengths.

#Now to plot the tree we use:
plot(NJtree, typ="unrooted", cex=0.7) # As our tree is unrooted I've changed type from "fan" to "unrooted"
title(expression("Neighbour-joining tree of "*italic(A.~fruticulosa)*" "))
#save the tree using:
write.tree((NJtree),file="NJ.tree_populations_alldata.tre")

# The correspondence between both analyses (pca and NJ) can be better assessed
# using colors based on PCs. We use this by using the colorplot function:
myCol <- colorplot(pca.1$scores, pca.1$scores, transp=TRUE, cex=2)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca.1$eig[1:28],2,1,2, posi="topright")

# Now the final pretty figure
plot(NJtree, typ="unrooted", show.tip=TRUE, cex=0.7)
tiplabels(pch=20, col=myCol, cex=6.5)
title(expression("Neighbour-joining tree of "*italic(A.~fruticulosa)*" "))

#Plot NJ tree with bootstrap values
library(gdsfmt)
library(SNPRelate) 
library(ggplot2)
library(ape)
setwd("C:/Users/matth/OneDrive/Documents/Bioinformatics")
fourdgvcf <- "4dg_noPA5.vcf" #read in VCF file
snpgdsVCF2GDS(fourdgvcf,"CEUexon2010_03genotype.gds",method ="biallelic.only")
genofileExon2010_03 <- snpgdsOpen("CEUexon2010_03genotype.gds")
set.seed(100)
ibs_Exon2010_03 <- snpgdsHCluster(snpgdsIBS(genofileExon2010_03,num.thread=2, autosome.only=FALSE))
rvExon2010_03 <- snpgdsCutTree(ibs_Exon2010_03)
treeExon2010_03 = rvExon2010_03$dendrogram
plot(rvExon2010_03$dendrogram,horiz=T, main ="CEU.exon.2010_03.genotypes.vcf SNP Tree" )
hcExon2010_03 = as.vector(rvExon2010_03$dendrogram)
install.packages("poppr")
library("poppr")
aboot(aa.D.ind, dist = nei.dist) #generates Bootstrap support values for tree