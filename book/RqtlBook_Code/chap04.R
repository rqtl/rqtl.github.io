######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 4: Single-QTL analysis
#######################################################################

######################################################################
# 4.1 Marker regression
######################################################################

# Load R/qtl and the hyper data
library(qtl)
data(hyper)

# Phenotype by genotype for two markers in the hyper data
par(mfrow=c(1,2))
plot.pxg(hyper, "D4Mit214")
plot.pxg(hyper, "D12Mit20")

# single-QTL scan by marker regression with the hyper data
out.mr <- scanone(hyper, method="mr")

# results for chromosome 12
out.mr[ out.mr$chr == 12, ]

# summary of out.mr
summary(out.mr, threshold=3)

# The biggest peak
max(out.mr)

# plot of marker regression results for chr 4 and 12
plot(out.mr, chr=c(4, 12), ylab="LOD score")

######################################################################
# 4.2 Interval mapping
######################################################################

###################################
# 4.2.1 Standard interval mapping
###################################

# Calculate genotype probabilities for the hyper data
hyper <- calc.genoprob(hyper, step=1, error.prob=0.001)

# illustrate the use of jittermap
hyper <- jittermap(hyper)

# Standard interval mapping for the hyper data
out.em <- scanone(hyper, method="em")

# Equivalent to the previous command
out.em <- scanone(hyper)

# Plot results of scanone with method="em"
plot(out.em, ylab="LOD score")

# Plot interval mapping and marker regression results together
plot(out.em, out.mr, chr=c(4, 12), col=c("blue", "red"),
     ylab="LOD score")

# Equivalent to the previous code
plot(out.em, chr=c(4, 12), col="blue", ylab="LOD score")
plot(out.mr, chr=c(4, 12), col="red", add=TRUE)

# A black/white version of the plot
plot(out.em, out.mr, chr=c(4, 12), col="black", lty=1:2,
     ylab="LOD score")

###################################
# 4.2.2 Haley--Knott regression
###################################

# Haley--Knott regression with the hyper data
out.hk <- scanone(hyper, method="hk")

# Plot results of EM and H-K together
plot(out.em, out.hk, chr=c(1,4,15), col=c("blue","red"),
     ylab="LOD score")

# Difference between H-K and EM results for hyper data
plot(out.hk - out.em, chr=c(1,4,15), ylim=c(-0.5, 1.0),
     ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

# Illustration of fancy math symbols in plots
plot(rnorm(100), rnorm(100), xlab=expression(hat(mu)[0]),
      ylab=expression(alpha^beta),
      main=expression(paste("Plot of ", alpha^beta, 
          " versus ", hat(mu)[0])))

###################################
# 4.2.3 Extended Haley--Knott regression
###################################

# extended Haley--Knott method with the hyper data
out.ehk <- scanone(hyper, method="ehk")

# Plot the three methods together
plot(out.em, out.hk, out.ehk, chr=c(1,4,15), ylab="LOD score",
     lty=c(1,1,2))

# Equivalent to the above
plot(out.em, chr=c(1,4,15), ylab="LOD score")
plot(out.hk, chr=c(1,4,15), col="blue", add=TRUE)
plot(out.ehk, chr=c(1,4,15), col="red", lty=2, add=TRUE)

# Plot differences between LOD scores
plot(out.hk - out.em, out.ehk - out.em, chr=c(1, 4, 15),
     col=c("blue", "red"), ylim=c(-0.5, 1),
     ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

###################################
# 4.2.4 Multiple imputation
###################################

# Perform multiple imputations with hyper data
hyper <- sim.geno(hyper, step=1, n.draws=64, error.prob=0.001)

# Interval mapping by multiple imputation with hyper data
out.imp <- scanone(hyper, method="imp")

# Plot imputation and em results for hyper data
plot(out.em, out.imp, chr=c(1,4,15), col=c("blue", "red"),
     ylab="LOD score")

# Difference between H-K and EM results for hyper data
plot(out.imp - out.em, chr=c(1,4,15), ylim=c(-0.5, 0.5),
     ylab=expression(LOD[IMP] - LOD[EM]))
abline(h=0, lty=3)

###################################
# 4.2.5 Comparison of methods
###################################

# Omit data from two chr 1 markers in listeria data
data(listeria)
mar2drop <- markernames(listeria, chr=1)[2:12]
listeria <- drop.markers(listeria, mar2drop)

# Run standard interval mapping 
listeria <- calc.genoprob(listeria, step=1, error.prob=0.001)
outl.em <- scanone(listeria, chr=1)
outl.hk <- scanone(listeria, chr=1, method="hk")
outl.ehk <- scanone(listeria, chr=1, method="ehk")
plot(outl.em, outl.hk, outl.ehk, ylab="LOD score", 
     lty=c(1,1,2))

# drop markers in hyper data having most missing data
data(hyper)
nt.mar <- ntyped(hyper, "mar")
mar2drop <- names(nt.mar[nt.mar < 92])
hyper.rev <- drop.markers(hyper, mar2drop)

# omit intermediate individuals' genotype data
nm.ind <- nmissing(hyper.rev)
ind2drop <- nm.ind > 0
for(i in 1:nchr(hyper))
  hyper.rev$geno[[i]]$data[ind2drop,] <- NA

# scanone with all individuals
hyper.rev <- calc.genoprob(hyper.rev, step=1, 
                           error.prob=0.001)
out1.em <- scanone(hyper.rev)
out1.hk <- scanone(hyper.rev, method="hk")
out1.ehk <- scanone(hyper.rev, method="ehk")

# plot LOD curves for EM, HK and eHK against one another
plot(out1.em, out1.hk, out1.ehk, ylab="LOD score", 
     lty=c(1,1,2))

# plot differences HK-EM and eHK-EM
plot(out1.hk-out1.em, out1.ehk-out1.em, col=c("blue", "red"),
     ylim=c(-0.1, 4), ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

# scanone with only the genotyped individuals
hyper.ex <- subset(hyper.rev, ind = !ind2drop)
out2.em <- scanone(hyper.ex)
out2.hk <- scanone(hyper.ex, method="hk")
out2.ehk <- scanone(hyper.ex, method="ehk")

# plot difference in LOD scores with fully genotyped individuals vs everyone
plot(out2.em-out1.em, out2.hk-out1.em, out2.ehk-out1.em,
     ylim=c(-0.1, 0.2), col=c("black", "green", "red"), 
     lty=c(1,3,2), ylab="Difference in LOD scores")
abline(h=0, lty=3)

# Time to run calc.genoprob
data(hyper)
system.time(hyper <- calc.genoprob(hyper, step=0.25, 
                                   error.prob=0.001))

# Time to run scanone
system.time(test.hk <- scanone(hyper, method="hk"))
system.time(test.ehk <- scanone(hyper, method="ehk"))
system.time(test.em <- scanone(hyper, method="em"))

# Time to do multiple imputation
system.time(hyper <- sim.geno(hyper, step=0.25, 
                              error.prob=0.001, n.draws=32))
system.time(test.imp <- scanone(hyper, method="imp"))

######################################################################
# 4.3 Significance thresholds
######################################################################

# Calculate genotype probabilities for the hyper data
data(hyper)
hyper <- calc.genoprob(hyper, step=1, error.prob=0.001)

# Code to illustrate a permutation test 
operm <- scanone(hyper, n.perm=1000, verbose=FALSE)

# First five permutation results 
operm[1:5]

# histogram of the permutation results
plot(operm)

# LOD thresholds
summary(operm, alpha=c(0.20, 0.05))

# define strata for hyper data
strat <- (ntyped(hyper) > 100)

# Code to illustrate a stratified permutation test 
operms <- scanone(hyper, n.perm=1000, perm.strata=strat,
                  verbose=FALSE)

summary(operms, alpha=c(0.20, 0.05))

# peaks meeting 10% significance level
summary(out.em, perms=operms, alpha=0.1)

# P-values for the peak LOD scores
summary(out.em, perms=operms, alpha=0.1, pvalues=TRUE)

# Upper confidence limit on true p-value.
binom.test(0, 1000)$conf.int

######################################################################
# 4.4 The X chromosome
######################################################################

###################################
# 4.4.1 Analysis
###################################

###################################
# 4.4.2 Significance thresholds
###################################

###################################
# 4.4.3 Example
###################################

# Access the iron data and make a plot
library(qtlbook)
data(iron)
plot(iron, pheno=1:2)

# scatterplot of log2(liver) and log2(spleen) in the iron data
plot(log2(liver) ~ log2(spleen), data=iron$pheno, 
     col=c("red", "blue")[iron$pheno$sex],
     pch=c(1,4)[iron$pheno$sex])

# Genome scan for the iron data
iron$pheno[,1] <- log2(iron$pheno[,1])
iron <- calc.genoprob(iron, step=1, error.prob=0.001)
out.liver <- scanone(iron)

# peaks with LOD>3
summary(out.liver, 3)

# autosome- and X-chr-specific permutation tests
operm.liver <- scanone(iron, n.perm=1000, perm.Xsp=TRUE,
                       verbose=FALSE)

# LOD thresholds
summary(operm.liver, alpha=0.05)

# Summary with p-values
summary(out.liver, perms=operm.liver, alpha=0.05, 
        pvalues=TRUE)

# Create strata
strat <- as.numeric(iron$pheno$sex) + iron$pheno$pgm*2
table(strat)

# Perform stratified permutation test
operm.liver.strat <- scanone(iron, n.perm=1000, perm.Xsp=TRUE,
                             perm.strata=strat, verbose=FALSE)

# LOD thresholds from stratified permutation test
summary(operm.liver.strat, alpha=0.05)

######################################################################
# 4.5 Interval estimates of QTL location
######################################################################

# LOD support and Bayes credible intervals for hyper, chr 4
lodint(out.em, 4, 1.5)
bayesint(out.em, 4, 0.95)

# Intervals for hyper, chr 4, expanded to the nearest flanking markers
lodint(out.em, 4, 1.5, expandtomarkers=TRUE)
bayesint(out.em, 4, 0.95, expandtomarkers=TRUE)

# bootstrap for QTL location
out.boot <- scanoneboot(hyper, chr=4, n.boot=1000)

# Get bootstrap confidence interval
summary(out.boot)

######################################################################
# 4.6 QTL effects
######################################################################

# find marker nearest 29.5 cM on chr 4
find.marker(hyper, 4, 29.5)

# effectplot with hyper data, marker D4Mit164
hyper <- sim.geno(hyper, n.draws=16, error.prob=0.001)
effectplot(hyper, mname1="D4Mit164")

# Get the estimates from effectplot
eff <- effectplot(hyper, mname1="D4Mit164", draw=FALSE)
eff

# Plot phenotypes against genotype at D4Mit164
plot.pxg(hyper, marker="D4Mit164")

######################################################################
# 4.7 Multiple phenotypes
######################################################################

# Add log(liver) and log(spleen) to the phenotypes
data(iron)
iron$pheno <- cbind(iron$pheno[,1:2], 
                    log2liver=log2(iron$pheno$liver),
                    log2spleen=log2(iron$pheno$spleen),
                    iron$pheno[,3:4])

# Genome scan with log(liver)
iron <- calc.genoprob(iron, step=1, error.prob=0.001)
out.logliver <- scanone(iron, pheno.col=3)

# Genome scan, referring to the phenotype by name.
out.logliver <- scanone(iron, pheno.col="log2liver")

# Genome scan, referring to the phenotype by name.
out.logliver <- scanone(iron, 
                        pheno.col=log2(iron$pheno$liver))

# Genome scan with all four phenotypes
out.all <- scanone(iron, pheno.col=1:4)

# Look at the beginning of out.all
out.all[1:5,]

# summary of out.all
summary(out.all, threshold=3, lodcolumn=4)

# summary of out.all with format="allpheno"
summary(out.all, threshold=3, format="allpheno")

# summary of out.all with format="allpeaks"
summary(out.all, threshold=3, format="allpeaks")

# plot of all LOD curves
plot(out.all, lodcolumn=1:2, col=c("blue", "red"), 
     chr=c(2, 7, 8, 9, 16), ylim=c(0,12.7), ylab="LOD score")
plot(out.all, lodcolumn=3:4, col=c("blue", "red"), lty=2,
     chr=c(2, 7, 8, 9, 16), add=TRUE)

# Permutations for all 4 iron phenotypes
operm.all <- scanone(iron, pheno.col=1:4, n.perm=1000, 
                     perm.Xsp=TRUE)

# LOD thresholds for all phenotypes
summary(operm.all, alpha=0.05)

# summary with 5% thresholds calculated from perms, and with pvalues
summary(out.all, format="allpeaks", perms=operm.all, 
        alpha=0.05, pvalues=TRUE)

# end of chap04.R
