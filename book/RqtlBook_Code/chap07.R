######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 7: Working with covariates
#######################################################################

######################################################################
# 7.1 Additive covariates
######################################################################

# Access the gutlength data
library(qtl)
library(qtlbook)
data(gutlength)

# Omit chromosome 15
gutlength <- subset(gutlength, chr = -15)

# Box plots of phenotype by sex and cross
boxplot(gutlength ~ sex*cross, data=gutlength$pheno, 
        horizontal=TRUE, xlab="Gut length (cm)",
        col=c("red","blue"))

# ANOVA relating gut length to sex and cross
anova(aov(gutlength ~ sex*cross, data=gutlength$pheno))

# Create indicators frommom and direction
cross <- as.numeric(pull.pheno(gutlength, "cross"))
frommom <- as.numeric(cross < 3)
forw <- as.numeric(cross == 1 | cross == 3)
gutlength$pheno$frommom <- frommom
gutlength$pheno$forw <- forw

# ANOVA relating gut length to sex and frommom, and forw
anova(aov(gutlength ~ sex*frommom*forw, data=gutlength$pheno))

# Make sex numeric 
sex <- as.numeric(pull.pheno(gutlength, "sex") == "M")

# Create a matrix for the cross covariate
crossX <- cbind(frommom, forw, frommom*forw)

# Paste sex and cross together
x <- cbind(sex, crossX)

# Genome scan with and without the covariates
gutlength <- calc.genoprob(gutlength, step=1, 
                           error.prob=0.001)
out.0 <- scanone(gutlength)
out.a <- scanone(gutlength, addcovar=x)

# Plot of the lod scores with and without the covariates
plot(out.0, out.a, col=c("blue", "red"), lty=1:2, 
     ylab="LOD score", alternate.chrid=TRUE)

# Plot of the lod score differences, with and without the covariates
plot(out.a - out.0, ylab="LOD w/ covar - LOD w/o covar", 
     ylim=c(-1, 1), alternate.chrid=TRUE)
abline(h=0, lty=2)

# Create strata regarding amount of genotyping
strat <- (nmissing(gutlength) < 50)

# Permutations with and without covariates, separately for autosomes and X chr
operm.0 <- scanone(gutlength, n.perm=1000, perm.Xsp=TRUE,
                   perm.strata=strat)
operm.a <- scanone(gutlength, addcovar=x, n.perm=1000, 
                   perm.Xsp=TRUE, perm.strata=strat)

# Combine results
out.both <- c(out.0, out.a, labels=c("nocovar", "covar"))
operm.both <- cbind(operm.0, operm.a, 
                    labels=c("nocovar", "covar"))

# LOD thresholds
summary(operm.both, 0.05)

# Summary of the main results
summary(out.both, perms=operm.both, format="allpeaks",
        alpha=0.2, pvalues=TRUE)

######################################################################
# 7.2 QTL $\times$ covariate interactions
######################################################################

# Genome scan with sex as an interactive covariate.
out.i <- scanone(gutlength, addcovar=x, intcovar=sex)

# plot of LODf and LODi
plot(out.i, out.i - out.a, ylab="LOD score", 
     col=c("blue", "red"), alternate.chrid=TRUE)

# Permutations with and without sex as interactive covariate
set.seed(54955149)
operm.a <- scanone(gutlength, addcovar=x, n.perm=1000, 
                   perm.Xsp=TRUE, perm.strata=strat)
set.seed(54955149)
operm.i <- scanone(gutlength, addcovar=x, intcovar=sex,
                   n.perm=1000, perm.Xsp=TRUE, 
                   perm.strata=strat)

# Combine the results
out.ia <- c(out.i, out.i - out.a, labels=c("f","i"))
operm.ia <- cbind(operm.i, operm.i - operm.a, 
                  labels=c("f","i"))

# Loci with large LOD_f 
summary(out.ia, perms=operm.ia, alpha=0.2, pvalues=TRUE)

# Loci with large LOD_i 
summary(out.ia, perms=operm.ia, alpha=0.2, pvalues=TRUE, 
        lodcolumn=2)

# Separate genome scan within each sex
out.m <- scanone(subset(gutlength, ind = sex==1), 
                 addcovar=crossX[sex==1,])
out.f <- scanone(subset(gutlength, ind = sex==0), 
                 addcovar=crossX[sex==0,])

# plot of LOD scores in males and females separately.
plot(out.m, out.f, col=c("blue", "red"), ylab="LOD score", 
     alternate.chrid=TRUE)

# Plot LOD(males) + LOD(females) - LOD_f
plot(out.m + out.f - out.i, ylim=c(-0.5,0.5), 
     ylab="LOD(males) + LOD(females) - LODf", 
     alternate.chrid=TRUE)
abline(h=0, lty=2)

# Permutations in each of males and females
operm.m <- scanone(subset(gutlength, ind = sex==1), 
                   addcovar=crossX[sex==1,], n.perm=1000,
                   perm.strata=strat[sex==1], perm.Xsp=TRUE)
operm.f <- scanone(subset(gutlength, ind = sex==0), 
                   addcovar=crossX[sex==0,], n.perm=1000,
                   perm.strata=strat[sex==0], perm.Xsp=TRUE)

# Combine the results
out.sexsp <- c(out.m, out.f, labels=c("male","female"))
operm.sexsp <- cbind(operm.m, operm.f, 
                     labels=c("male","female"))

# LOD thresholds
summary(operm.sexsp, 0.05)

# LOD thresholds
summary(out.sexsp, perms=operm.sexsp, alpha=0.2, 
        pvalues=TRUE, format="allpeaks")

# Effect plots for four selected chromosomes
gutlength <- sim.geno(gutlength, n.draws=128, step=1,
                      error.prob=0.001)
par(mfrow=c(2,2))
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="4@65", main="Chromosome 4", 
           add.legend=FALSE)
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="5@22", main="Chromosome 5", 
           add.legend=FALSE)
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="18@48.2", main="Chromosome 18", 
           add.legend=FALSE)
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="X@58", main="X chromosome", 
           add.legend=FALSE)

######################################################################
# 7.3 Covariates with non-normal phenotypes
######################################################################

# load the nf1 data set
data(nf1)
nf1 <- drop.nullmarkers(nf1)

# Proportion affected
mean(pull.pheno(nf1,"affected"))
tapply(pull.pheno(nf1,"affected"), 
       pull.pheno(nf1,"from.mom"), mean)

# chi-square test
chisq.test(pull.pheno(nf1,"affected"), 
           pull.pheno(nf1,"from.mom"))

# genome scans for nf1
nf1 <- calc.genoprob(nf1, step=1, error.prob=0.001)
from.mom <- pull.pheno(nf1,"from.mom")
out.a <- scanone(nf1, model="binary", addcovar=from.mom)
out.i <- scanone(nf1, model="binary", addcovar=from.mom,
                 intcovar=from.mom)

# Permutation tests
set.seed(1310709)
operm.a <- scanone(nf1, model="binary", addcovar=from.mom,
                   n.perm=1000)
set.seed(1310709)
operm.i <- scanone(nf1, model="binary", addcovar=from.mom,
                   intcovar=from.mom, n.perm=1000)

# Combine the results
out.all <- c(out.i, out.a, out.i-out.a, labels=c("f","a","i"))
operm.all <- cbind(operm.i, operm.a, operm.i - operm.a,
                   labels=c("f","a","i"))

# Plot nf1 results
plot(out.all, lod=1:3, ylab="LOD score")

summary(out.all, perms=operm.all, alpha=0.2, pvalues=TRUE)

pchisq(0.703 * 2 * log(10), 1, lower=FALSE)

# nf1 scans split by parent-of-origin
out.frommom <- scanone(subset(nf1, ind=(from.mom==1)),
                       model="binary")
out.fromdad <- scanone(subset(nf1, ind=(from.mom==0)),
                       model="binary")

# Permutation tests
operm.frommom <- scanone(subset(nf1, ind=(from.mom==1)),
                       model="binary", n.perm=1000)
operm.fromdad <- scanone(subset(nf1, ind=(from.mom==0)),
                       model="binary", n.perm=1000)

# Combine the results
out.bypoo <- c(out.frommom, out.fromdad, 
               labels=c("mom","dad"))
operm.bypoo <- cbind(operm.frommom, operm.fromdad,
                     labels=c("mom","dad"))

# Plot nf1, split by parent-of-origin, results
plot(out.bypoo, lod=1:2, col=c("red", "blue"), 
     ylab="LOD score")

# Summary of the results
summary(out.bypoo, perms=operm.bypoo, alpha=0.2, pvalues=TRUE,
        format="allpeaks")

# Plot proportion affected as a function of two-locus genotypes
nf1 <- sim.geno(nf1, n.draws=128, step=1, error.prob=0.001)
par(mfrow=c(1,2))
effectplot(nf1, mname1="NPcis", mark1=1-from.mom, 
           geno1=c("Mom", "Dad"), mname2="15@13", 
           ylim=c(0,1))
effectplot(nf1, mname1="NPcis", mark1=1-from.mom, 
           geno1=c("Mom", "Dad"), mname2="19@0", ylim=c(0,1))

######################################################################
# 7.4 Composite interval mapping
######################################################################

# Reload hyper and rerun genome scan
data(hyper)
hyper <- calc.genoprob(hyper, step=1, error.prob=0.001)
out <- scanone(hyper)

# Genotype data for marker near peak LOD score
mar <- find.marker(hyper, 4, 29.5)
g <- pull.geno(hyper)[,mar]
sum(is.na(g))

# Imputate genotype data
g <- pull.geno(fill.geno(hyper))[,mar]
sum(is.na(g))

# How to recode the genotype data for an intercross [don't run here]
g <- cbind(as.numeric(g==1), as.numeric(g==2))

# genome scan with the marker as an additive covariate
out.ag <- scanone(hyper, addcovar=g)

# Plot the results
plot(out, out.ag, col=c("blue", "red"), ylab="LOD score")

# Permutation test
strat <- (ntyped(hyper) > 100)
operm.ag <- scanone(hyper, addcovar=g, chr=-4, 
                    perm.strata=strat, n.perm=1000)

# LOD thresholds
summary(operm.ag, alpha=c(0.2, 0.05))

# Summary of results
summary(out.ag, perms=operm.ag, alpha=0.2, pvalues=TRUE)

# genome scan with the marker as an interactive covariate
out.ig <- scanone(hyper, addcovar=g, intcovar=g)

# Plot the interaction LOD scores
plot(out.ig - out.ag, ylab="interaction LOD score")

# run CIM with different window sizes
out.cim.20 <- cim(hyper, n.marcovar=3, window=20)
out.cim.40 <- cim(hyper, n.marcovar=3, window=40)
out.cim.inf <- cim(hyper, n.marcovar=3, window=Inf)

# Plot the CIM results
chr <- c(1, 2, 4, 6, 15)
par(mfrow=c(3,1))
plot(out, out.cim.20, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = 20 cM")
add.cim.covar(out.cim.20, chr=chr, col="green")
plot(out, out.cim.40, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = 40 cM")
add.cim.covar(out.cim.40, chr=chr, col="green")
plot(out, out.cim.inf, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = Inf")
add.cim.covar(out.cim.inf, chr=chr, col="green")

# end of chap07.R
