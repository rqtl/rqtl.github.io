######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap03.R"
#
# Chapter 3: Data checking
#######################################################################

######################################################################
# Figure 3.3
#
#  Plots of the average phenotype against index (left panel) and
#  against a randomized index (right panel) for the ch3a
#  data.
######################################################################
set.seed(55767947)
par(mar=c(4.1,4.1,0.1,1.1))
par(mfrow=c(1,2), las=1, cex=0.8)
means <- apply(ch3a$pheno, 1, mean)
plot(means)
plot(sample(means), xlab="Random index", ylab="means")



######################################################################
# Figure 3.4
#
#  Histogram of the proportion of markers with identical
#  genotypes for each pair of individuals in the ch3a
#  data.
######################################################################
data(ch3a)
cg <- comparegeno(ch3a)
par(mar=c(4.1,1.1,0.1,1.1))
hist(cg, breaks=seq(0, 1, len=201), main="", yaxt="n", ylab="",
     xlab="Proportion of identical genotypes")
rug(cg)



######################################################################
# Figure 3.7
#
#  Genetic map, as estimated from the ch3c
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
nm <- est.map(ch3c, error.prob=0.001)
plot(nm, main="")



######################################################################
# Figure 3.11
#
#  The genetic map in the ch3c data (after considerable
#  revisions in marker order) plotted against the map estimated from
#  the data.  For each chromosome, the line on the left is the map in
#  the data, and the line on the right is the
#  map estimated from the data; line segments connect the positions for each
#  marker. 
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
nm <- est.map(ch3c, error.prob=0.001, verbose=FALSE)
plot.map(ch3c, nm)



######################################################################
# Figure 3.12
#
#  Chromosome 16 genotypes for selected individuals in the
#  hyper data.  Open and closed circles are homozygous and
#  heterozygous genotypes, respectively. Possible genotyping errors are
#  flagged with red squares; inferred crossovers are indicated with
#  blue x's.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1),las=1)
plot.geno(hyper, 16, top$id[top$chr==16], cutoff=5,
          min.sep=4, main="")



######################################################################
# Figure 3.13
#
#  The observed number of crossovers for each individual in the
#  hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1),las=1)
nxo <- countXO(hyper)
plot(nxo, ylab="No. crossovers")



######################################################################
# Figure 3.14
#
#  The proportion of missing genotype information in the
#  hyper data.  Results by the entropy and variance versions are
#  shown in blue and red, respectively.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot.info(hyper, col=c("blue", "red"), ylab="Missing information", 
          main="")



######################################################################
# Figure 3.15
#
#  Histogram of the number of individuals with missing genotypes
#  at the markers in the hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1),las=1)
hist(nmissing(hyper, what="mar"), breaks=50,
     xlab="No. missing genotypes", main="")


# end of fig03.R
