######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 11: Case study II
#######################################################################

######################################################################
# 11.1 Diagnostics
######################################################################

# Load the data
library(qtl)
library(qtlbook)
data(trout)

# cross type
class(trout)

# Example of reading DH data (not run here)
trout <- read.cross("csv", file="trout.csv", 
                    genotypes=c("C","O"), alleles=c("C","O"))

# Example of how to change the cross type (not run here)
class(trout)[1] <- "dh"

# summary of the data
summary(trout)

# summary plot of the data
plot(trout)

# boxplot of tth vs egg source
boxplot(tth ~ female, data=trout$pheno, 
        xlab="Female", ylab="Time to hatch")

# anova to compare tth across females
anova(aov(tth ~ female, data=trout$pheno))

# markers with aberrant segregation patterns 
gt <- geno.table(trout)
gt[gt$P.value < 0.05/totmar(trout),]

# estimated pairwise recombination fractions
trout <- est.rf(trout)
plot.rf(trout, alternate.chrid=TRUE)

# plot rf for selected pairs of linkage groups
plot.rf(trout, chr=c(2,29, 10,18, 12,16, 14,20, 27,31),
        alternate.chrid=TRUE)

# table of two-locus genotypes
mar <- find.marker(trout, c(10,18), c(0,0))
geno.crosstab(trout, mar[1], mar[2])

# swap genotypes on the second linkage group of each pair
trout$geno[["29"]]$data <- 3 - trout$geno[["29"]]$data
trout$geno[["18"]]$data <- 3 - trout$geno[["18"]]$data
trout$geno[["16"]]$data <- 3 - trout$geno[["16"]]$data
trout$geno[["20"]]$data <- 3 - trout$geno[["20"]]$data
trout$geno[["31"]]$data <- 3 - trout$geno[["31"]]$data

# move markers from lg 18 to lg 10
lg18mar <- markernames(trout, 18)
for(i in lg18mar)
  trout <- movemarker(trout, i, 10)

# Change the name of the merged linkage group
nam <- names(trout$geno)
nam[nam=="10"] <- "10.18"
names(trout$geno) <- nam

# try all possible marker orders on lg 10.18
rip <- ripple(trout, chr="10.18", window=6, method="lik", 
              error.prob=0.01, map.function="kosambi",
              verbose=FALSE)

# summary of ripple results
summary(rip)

# switch order of markers on lg 10.18
trout <- switch.order(trout, "10.18", rip[2,], 
                      error.prob=0.01, map.function="kosambi")

# look at map of lg 10.18
pull.map(trout, chr="10.18")

# fix the other pairs of linkage groups
lg29mar <- markernames(trout, 29)
for(i in lg29mar)
  trout <- movemarker(trout, i, 2)
lg16mar <- markernames(trout, 16)
for(i in lg16mar)
  trout <- movemarker(trout, i, 12)
trout <- switch.order(trout, chr=12, c(1:8,10,11,9), 
                      error.prob=0.01, map.function="kosambi")
lg20mar <- markernames(trout, 20)
for(i in lg20mar)
  trout <- movemarker(trout, i, 14)
trout <- switch.order(trout, chr=14, error.prob=0.01,
                      c(10:1,14,16,15,17,18,13:11),
                      map.function="kosambi")
lg31mar <- markernames(trout, 31)
for(i in lg31mar)
  trout <- movemarker(trout, i, 27)
trout <- switch.order(trout, chr=27, error.prob=0.01,
                      c(1:6,8,7,9:13,15:18,14,19),
                      map.function="kosambi")

# Change the names of the merged linkage groups
nam <- names(trout$geno)
nam[nam=="2"] <- "2.29"
nam[nam=="12"] <- "12.16"
nam[nam=="14"] <- "14.20"
nam[nam=="27"] <- "27.31"
names(trout$geno) <- nam

# estimated r.f. for selected linkage groups
plot.rf(trout, chr=c(13, 19, "12.16", 25))

# reestimate map and plot it against that in data
newmap <- est.map(trout, error.prob=0.01, verbose=FALSE,
                  map.function="kosambi")
plot.map(trout, newmap, alternate.chrid=TRUE)

# Replace map with the newly estimated one
trout <- replace.map(trout, newmap)

######################################################################
# 11.2 Initial QTL analyses
######################################################################

# Calculate conditional QTL genotype probabilites
trout <- calc.genoprob(trout, step=1, err=0.01, 
                       map.function="kosambi")

# set up covariates
female <- pull.pheno(trout, "female")
lev <- levels(female)
nlev <- length(lev)
femcov <- matrix(0, nrow=nind(trout), ncol=nlev-1)
colnames(femcov) <- lev[-1]
for(i in 2:nlev)
  femcov[female==lev[i],i-1] <- 1

# genome scan with MCE as additive covariates
out <- scanone(trout, method="hk", addcovar=femcov)

# plot of lod curves
plot(out, ylab="LOD score", alternate.chrid=TRUE)

# permutation test
set.seed(523938)
operm <- scanone(trout, method="hk", addcovar=femcov, 
                 n.perm=4000)

# summary of genome scan
summary(out, perms=operm, alpha=0.1, pvalues=TRUE)

# genome scan controlling for QTL on lg 8
qtl <- makeqtl(trout, 8, 9.11, what="prob")
out.c8 <- addqtl(trout, qtl=qtl, method="hk", covar=femcov)

# summary of scan controlling for the QTL on lg 8
summary(out.c8, perms=operm, alpha=0.1, pvalues=TRUE)

# plot scan w/ and w/o controlling for lg 8 locus
plot(out, out.c8, chr = -8, col=c("blue","red"), 
     ylab="LOD score", alternate.chrid=TRUE)

# combined model; refinement of positions
qtl <- makeqtl(trout, c(8, 9, "10.18", "14.20", 17, 24),
               c(9.11, 21, 10, 0, 20, 44), what="prob")
rqtl <- refineqtl(trout, qtl=qtl, covar=femcov, 
                  method="hk", verbose=FALSE)

# print results
options(width=64)
rqtl

# drop-one-QTL analysis 
summary(fitqtl(trout, qtl=rqtl, covar=femcov, method="hk"),
        pvalues=FALSE)

# test for interactions
addint(trout, qtl=rqtl, qtl.only=TRUE, method="hk", 
       covar=femcov, pvalues=FALSE)

# 2d 2-qtl permutation test
operm2 <- scantwo(trout, method="hk", addcovar=femcov, 
                  n.perm=1000)

# significance thresholds from 2d scan
summary(operm2)

# perform 2d, 2-qtl scan
out2 <- scantwo(trout, method="hk", addcovar=femcov,
                incl.markers=TRUE)

# summary of scantwo output
summary(out2, perms=operm2, alpha=0.05)

# 2d-scan plots
plot(out2, lower="cond-int", chr=8)
plot(out2, lower="cond-int", chr=c(7, "12.16", 13))

# find the markers near the inferred QTL pairs
mar7 <- find.marker(trout, 7, c(11.52, 16))
find.markerpos(trout, mar7)
mar8 <- find.marker(trout, 8, c(9.11, 29))
find.markerpos(trout, mar8)
mar12.16 <- find.marker(trout, "12.16", c(2, 6.95))
find.markerpos(trout, mar12.16)
mar13 <- find.marker(trout, 13, c(0, 23))
find.markerpos(trout, mar13)

# plot phenotype against genotype
par(mfrow=c(2,2))
plot.pxg(trout, marker=mar7)
plot.pxg(trout, marker=mar8)
plot.pxg(trout, marker=mar12.16)
plot.pxg(trout, marker=mar13)

# bring all QTL together and fit multiple-QTL model
qtl <- makeqtl(trout, c(7,7,8,8,"10.18","12.16","12.16",13,13,24),
               c(11.52,16,9.11,29,11,2,6.95,0,23,46), what="prob")
summary(fitqtl(trout, qtl=qtl, covar=femcov, method="hk", 
               formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                         Q1+Q2+Q3*Q4+Q5+Q6*Q7+Q8*Q9+Q10), 
        pvalues=FALSE)

# drop QTL on 12.16 and 13 and refine QTL positions
qtl2 <- dropfromqtl(qtl, 6:9)
qtl2 <- refineqtl(trout, qtl=qtl2, covar=femcov, method="hk", 
                  formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                  Q1+Q2+Q3*Q4+Q5+Q6, verbose=FALSE)

# drop-one-QTL analysis
summary(fitqtl(trout, qtl=qtl2, covar=femcov, method="hk", 
               formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                         Q1+Q2+Q3*Q4+Q5+Q6),
        pvalues=FALSE)

# drop 2nd QTL on 7 and refine QTL positions
qtl3 <- dropfromqtl(qtl2, 2)
qtl3 <- refineqtl(trout, qtl=qtl3, covar=femcov, method="hk", 
                  formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                  Q1+Q2*Q3+Q4+Q5, verbose=FALSE)

# drop-one-QTL analysis
summary(fitqtl(trout, qtl=qtl3, covar=femcov, method="hk", 
               formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                         Q1+Q2*Q3+Q4+Q5),
        pvalues=FALSE)

# Omit locus on lg 7 and rerun refineqtl
qtl4 <- dropfromqtl(qtl3, 1)
qtl4 <- refineqtl(trout, qtl=qtl4, covar=femcov, method="hk", 
                  formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                  Q1*Q2+Q3+Q4, verbose=FALSE)

# drop-one-QTL analysis
summary(fitqtl(trout, qtl=qtl4, covar=femcov, method="hk", 
               formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                         Q1*Q2+Q3+Q4),
        pvalues=FALSE)

# calculate penalties 
print(pen <- calc.penalties(operm2))

# stepwise analysis with no interactions
stepout1 <- stepwiseqtl(trout, covar=femcov, penalties=pen,
                        method="hk", additive.only=TRUE, 
                        verbose=FALSE)

# print result
stepout1

# drop-one analysis
summary(fitqtl(trout, qtl=stepout1, covar=femcov, 
               method="hk"), pvalues=FALSE)

# stepwise analysis allowing interactions
stepout2 <- stepwiseqtl(trout, covar=femcov, penalties=pen,
                       method="hk", verbose=FALSE)

# print result
stepout2

# estimated effects
summary(fitqtl(trout, qtl=stepout2, covar=femcov, 
               method="hk", dropone=FALSE, get.ests=TRUE,
               formula=y~TL3+SP1+SP2+SP3+SP4+SP5+SP6+
                         Q1*Q2+Q3+Q4))

######################################################################
# 11.3 QTL $\times$ covariate interactions
######################################################################

# genome scan with QTL x MCE interactions
outi <- scanone(trout, method="hk", addcovar=femcov,
                intcovar=femcov)

# combine the three sets of LOD curves
outi <- c(out, outi, outi-out, labels=c("a","f","i"))

# plot LOD scores
plot(outi, lod=1:3, ylab="LOD score", alternate.chrid=TRUE)

# permutation test with MCE as interactive covariate
set.seed(523938)
opermi <- scanone(trout, method="hk", addcovar=femcov, 
                  intcovar=femcov, n.perm=4000)

# Combine the permutation test results
opermi <- cbind(operm, opermi, opermi-operm, 
                labels=c("a","f","i"))

# thresholds
summary(opermi)

# pointwise threshold for QTL x MCE interaction
qchisq(0.95, 7) / (2 * log(10))

# summary of scanone results 
summary(outi, perms=opermi, alpha=0.05, pvalues=TRUE,
        format="allpheno")

# make QTL object with interacting loci on linkage group 8
qtl <- makeqtl(trout, c("8", "8"), c(9.1, 28), what="prob")

# create formulas
addform <- paste("y~Q1*Q2+Q3+", 
                 paste(colnames(femcov), collapse="+"),
                 sep="")
addform

intform <- paste("y~Q1*Q2+Q3+", 
                 paste("Q3", colnames(femcov),
                       sep="*", collapse="+"),
                 sep="")
intform

# adjust for lg 8 loci
out.aq <- addqtl(trout, qtl=qtl, method="hk", covar=femcov, 
                 formula=addform)
outi.aq <- addqtl(trout, qtl=qtl, method="hk", covar=femcov,
                  formula=intform)

# combine the three sets of LOD scores
outi.aq <- c(out.aq, outi.aq, outi.aq - out.aq, 
             labels=c("a","f","i"))

# summary of results
summary(outi.aq, perms=opermi, alpha=0.05, pvalues=TRUE,
        format="allpheno")

# plot LOD scores
plot(outi.aq, lod=1:3, ylab="LOD score", alternate.chrid=TRUE,
     chr=c(8, 9, "10.18", "12.16", 17, 24))

# add qtl to object
qtl <- addtoqtl(trout, qtl, c(9, "10.18", 17, 24), 
                c(30, 28.7, 18, 44))

# create full model formula
fullform <- paste("y~Q1*Q2+Q3+Q4+Q5+Q6", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q3", colnames(femcov), sep=":", 
                        collapse="+"),
                  paste("Q4", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")

# formulas with one QTL omitted
form.m9 <- paste("y~Q1*Q2+Q4+Q5+Q6", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q4", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")
form.m1018 <- paste("y~Q1*Q2+Q3+Q5+Q6", 
                    paste(colnames(femcov), collapse="+"),
                    paste("Q3", colnames(femcov), sep=":",
                          collapse="+"), sep="+")
form.m17 <- paste("y~Q1*Q2+Q3+Q4+Q6", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q3", colnames(femcov), sep=":", 
                        collapse="+"), 
                  paste("Q4", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")
form.m24 <- paste("y~Q1*Q2+Q3+Q4+Q5", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q3", colnames(femcov), sep=":", 
                        collapse="+"), 
                  paste("Q4", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")

# formulas with QTL x MCE interaction omitted
form.m9int <- paste("y~Q1*Q2+Q3+Q4+Q5+Q6", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q4", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")
form.m1018int <- paste("y~Q1*Q2+Q3+Q4+Q5+Q6", 
                       paste(colnames(femcov), collapse="+"),
                       paste("Q3", colnames(femcov), sep=":",
                             collapse="+"), sep="+")

# refine qtl positions
qtl <- refineqtl(trout, qtl=qtl, formula=fullform, 
                 method="hk", covar=femcov, verbose=FALSE)

# fit full model
full <- fitqtl(trout, qtl=qtl, formula=fullform, method="hk",
               covar=femcov, dropone=FALSE)

# summary of full model
summary(full, pvalues=FALSE)

# pull out lod score for full model
print(fulllod <- full$lod)

# fit all other models
m9 <- fitqtl(trout, qtl=qtl, formula=form.m9, method="hk",
             covar=femcov, dropone=FALSE)
m1018 <- fitqtl(trout, qtl=qtl, formula=form.m1018, 
                method="hk", covar=femcov, dropone=FALSE)
m17 <- fitqtl(trout, qtl=qtl, formula=form.m17, method="hk",
              covar=femcov, dropone=FALSE)
m24 <- fitqtl(trout, qtl=qtl, formula=form.m24, method="hk",
              covar=femcov, dropone=FALSE)
m9int <- fitqtl(trout, qtl=qtl, formula=form.m9int, 
                method="hk", covar=femcov, dropone=FALSE)
m1018int <- fitqtl(trout, qtl=qtl, formula=form.m1018int, 
                   method="hk", covar=femcov, dropone=FALSE)

# LOD scores for omitting loci on lg 17 and 24
fulllod - m17$lod
fulllod - m24$lod

# full lod for lg 9 and 10.18
fulllod - m9$lod
fulllod - m1018$lod

# interation LOD scores
fulllod - m9int$lod
fulllod - m1018int$lod

# refineqtl for all of the submodels
qtl.m9 <- refineqtl(trout, qtl=qtl, formula=form.m9,
                    method="hk", covar=femcov,
                    verbose=FALSE)
qtl.m1018 <- refineqtl(trout, qtl=qtl, formula=form.m1018,
                       method="hk", covar=femcov,
                       verbose=FALSE)
qtl.m17 <- refineqtl(trout, qtl=qtl, formula=form.m17,
                     method="hk", covar=femcov,
                     verbose=FALSE)
qtl.m24 <- refineqtl(trout, qtl=qtl, formula=form.m24,
                     method="hk", covar=femcov,
                     verbose=FALSE)
qtl.m9int <- refineqtl(trout, qtl=qtl, formula=form.m9int,
                       method="hk", covar=femcov,
                       verbose=FALSE)
qtl.m1018int <- refineqtl(trout, qtl=qtl, method="hk",
                          formula=form.m1018int,covar=femcov,
                          verbose=FALSE)

# fit reduced models with QTL in refined positions
m9r <- fitqtl(trout, qtl=qtl.m9, formula=form.m9, 
              method="hk", covar=femcov, dropone=FALSE)
m1018r <- fitqtl(trout, qtl=qtl.m1018, formula=form.m1018,
                 method="hk", covar=femcov, dropone=FALSE)
m17r <- fitqtl(trout, qtl=qtl.m17, formula=form.m17, 
               method="hk", covar=femcov, dropone=FALSE)
m24r <- fitqtl(trout, qtl=qtl.m24, formula=form.m24, 
               method="hk", covar=femcov, dropone=FALSE)
m9intr <- fitqtl(trout, qtl=qtl.m9int, formula=form.m9int,
                 method="hk", covar=femcov, dropone=FALSE)
m1018intr <- fitqtl(trout, qtl=qtl.m1018int, method="hk",
                    formula=form.m1018int, covar=femcov, 
                    dropone=FALSE)

# revised LOD scores for omitting loci on lg 17 and 24
fulllod - m17r$lod
fulllod - m24r$lod

# revised full lod for lg 9 and 10.18
fulllod - m9r$lod
fulllod - m1018r$lod
fulllod - m9intr$lod
fulllod - m1018intr$lod

# plot LOD profiles for our 6-QTL model
plotLodProfile(qtl, col=c("blue","red",rep("black",4)),
               ylab="Profile LOD score")

# CI for QTL on linkage group 10.18
lodint(qtl, qtl.index=4)
bayesint(qtl, qtl.index=4)

# genetic map of markers on linkage group 10.18
pull.map(trout, "10.18")

# drop linkage group 8 QTL
qtl.m8 <- dropfromqtl(qtl, 1:2)

# form model formula
theformula <- paste("y~Q1+Q2+Q3+Q4+Q5*Q6", 
                  paste(colnames(femcov), collapse="+"),
                  paste("Q1", colnames(femcov), sep=":", 
                        collapse="+"),
                  paste("Q2", colnames(femcov), sep=":", 
                        collapse="+"), sep="+")

# 2d scan on linkage group 8
out.ap <- addpair(trout, chr="8", qtl=qtl.m8, covar=femcov,
                  formula=theformula, method="hk", 
                  incl.markers=TRUE, verbose=FALSE)

# Plot 2d profile with 1.5-LOD contour
plot(out.ap, contours=1.5)

# estimated effects
summary(fitqtl(trout, qtl=qtl, formula=fullform, covar=femcov,
               method="hk", dropone=FALSE, get.ests=TRUE))

######################################################################
# 11.4 Discussion
######################################################################

# end of chap11.R
