######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 3: Data checking
#######################################################################

######################################################################
# 3.1 Phenotypes
######################################################################

# Load the qtl and qtlbook packages and the ch3a data
library(qtl)
library(qtlbook)
data(ch3a)

# Code to make histograms of each phenotype in the ch3a data
par(mfrow=c(3,2))
for(i in 1:5)
  plot.pheno(ch3a, pheno.col=i)

# Code to create a scatterplot matrix of the phenotypes
# in the ch3a data
pairs(jitter( as.matrix(ch3a$pheno) ), cex=0.6, las=1)

# Remove 0 phenotypes
ch3a$pheno[ch3a$pheno == 0] <- NA

# Code to create a plots of each individual's average phenotype
# against its index and then against randomized indices
par(mfrow=c(1,2), las=1, cex=0.8)
means <- apply(ch3a$pheno, 1, mean)
plot(means)
plot(sample(means), xlab="Random index", ylab="means")

######################################################################
# 3.2 Segregation distortion
######################################################################

# Segregation distortion in the ch3b data
data(ch3b)
gt <- geno.table(ch3b)
gt[ gt$P.value < 1e-7, ]

# Segregation distortion in the listeria data
data(listeria)
gt <- geno.table(listeria)
p <- gt$P.value
gt[ !is.na(p) & p < 0.01, ]

######################################################################
# 3.3 Compare individuals' genotypes
######################################################################

# Histogram of the proportion of markers with identical 
# genotypes for each pair of individuals in the ch3a data.
data(ch3a)
cg <- comparegeno(ch3a)
hist(cg, breaks=200, 
     xlab="Proportion of identical genotypes")
rug(cg)

# Identify the relevant pairs
which(cg > 0.9, arr.ind=TRUE)

######################################################################
# 3.4 Check marker order
######################################################################

###################################
# 3.4.1 Pairwise recombination fractions
###################################

# Calculate recombination fractions between all pairs of 
# markers in the ch3c data
data(ch3c)
ch3c <- est.rf(ch3c)

# Get details on the markers with alleles potentially swapped
checkAlleles(ch3c)

# Look more carefully at the genotypes on chr 1,
# in order to figure out which of the markers
# has problems
pull.map(ch3c, 1)
geno.crosstab(ch3c, "c1m3", "c1m4")
geno.crosstab(ch3c, "c1m3", "c1m5")
geno.crosstab(ch3c, "c1m4", "c1m5")

# Swap the alleles at marker 2 and put the
# data back in the ch3c object.
g <- pull.geno(ch3c, 1)
g[,"c1m3"] <- 4 - g[,"c1m3"]
ch3c$geno[[1]]$data <- g

# Swap the alleles at marker 2 on chr 7
g <- pull.geno(ch3c, chr=7)
g[,"c7m2"] <- 4 - g[,"c7m2"]
ch3c$geno[[7]]$data <- g

# Rerun est.rf and checkAlleles
ch3c <- est.rf(ch3c)
checkAlleles(ch3c)

# Plot the recombination fractions for the ch3c data
plot.rf(ch3c, alternate.chrid=TRUE)

# Plot rec fracs for selected chromosomes
plot.rf(ch3c, chr=c(1, 7, 12, 13, 15))

# Estimate genetic map and plot it
nm <- est.map(ch3c, error.prob=0.001)
plot(nm)

# Get names of problem markers and then move them to the
# chromosomes to which they seem to belong
ch3c <- movemarker(ch3c, find.marker(ch3c, 7,index=4), 15)
ch3c <- movemarker(ch3c, find.marker(ch3c,12,index=2),  1)
ch3c <- movemarker(ch3c, find.marker(ch3c,12,index=1),  7)
ch3c <- movemarker(ch3c, find.marker(ch3c,13,index=1), 12)
ch3c <- movemarker(ch3c, find.marker(ch3c,15,index=5), 12)

# Plot the recombination fractions for the hyper data.
data(hyper)
hyper <- est.rf(hyper)
plot.rf(hyper, alternate.chrid=TRUE)

###################################
# 3.4.2 Rippling marker order
###################################

# Assess marker orders on chr 1 in the ch3c data
rip <- ripple(ch3c, 1, 5)

# Summary of the ripple results
summary(rip)

# Switch order of markers on chr 1
ch3c <- switch.order(ch3c, 1, rip[2,])

# ripple via maximum likelihood
rip <- ripple(ch3c, 1, 3, method="likelihood",
              error.prob=0.001, verbose=FALSE)
summary(rip)

# Apply ripple to all chromosomes
rip <- vector("list", nchr(ch3c))
names(rip) <- names(ch3c$geno)
for(i in names(ch3c$geno))
  rip[[i]] <- ripple(ch3c, i, 7, verbose=FALSE)

# Get improvement in number of obligate crossovers
dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])
dif.nxo

# Switch orders on chromosomes showing possible improvement
for(i in names(ch3c$geno)) {
  if(dif.nxo[i] > 0)
    ch3c <- switch.order(ch3c, i, rip[[i]][2,])
}

# Repeat ripple with method="countxo"
for(i in names(ch3c$geno))
  rip[[i]] <- ripple(ch3c, i, 7, verbose=FALSE)
dif.nxo <- sapply(rip, function(a) a[1,ncol(a)]-a[2,ncol(a)])

any(dif.nxo > 0)

# Ripple all chromosomes by maximum likelihood
for(i in names(ch3c$geno))
  rip[[i]] <- ripple(ch3c, i, 3, method="likelihood", 
                     error.prob=0.001, verbose=FALSE)
lod <- sapply(rip, function(a) a[2, ncol(a)-1])

lod[lod > 0]

# summary of ripple results for X chromosome
summary(rip[["X"]])

# rerun est.rf and plot.rf
ch3c <- est.rf(ch3c)
plot.rf(ch3c, chr=c(1, 7, 12, 13, 15))

###################################
# 3.4.3 Estimate genetic map
###################################

# estimate map for ch3c data and plot it against that in the data
nm <- est.map(ch3c, error.prob=0.001, verbose=FALSE)
plot.map(ch3c, nm)

# replace map
ch3c <- replace.map(ch3c, nm)

######################################################################
# 3.5 Identifying genotyping errors
######################################################################

# calculate genotyping error LOD scores
data(hyper)
newmap <- est.map(hyper, error.prob=0.01)
hyper <- replace.map(hyper, newmap)
hyper <- calc.errorlod(hyper)

# The top error LOD scores for the hyper data
top <- top.errorlod(hyper, cutoff=5)
top

# plot genotypes for hyper data, chr 16
plot.geno(hyper, 16, top$id[top$chr==16], cutoff=5)

######################################################################
# 3.6 Counting crossovers
######################################################################

nxo <- countXO(hyper)
plot(nxo, ylab="No. crossovers")

# Individuals with many crossovers
nxo[nxo>25]

# ave no. XOs for first 92 and last 158 individuals
mean(nxo[1:92])
mean(nxo[-(1:92)])

# No. XOs by chromosome for individual 56.
countXO(hyper, bychr=TRUE)[56,]

######################################################################
# 3.7 Missing genotype information
######################################################################

# plot of missing genotype information in hyper data
plot.info(hyper, col=c("blue", "red"))

# save plot info results in an object
z <- plot.info(hyper, step=0)
z[ z[,1]==14, ]

# Plot of Number of missing genotypes at each marker
hist(nmissing(hyper, what="mar"), breaks=50)

# end of chap03.R
