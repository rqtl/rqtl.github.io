######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 10: Case study I
#######################################################################

######################################################################
# 10.1 Diagnostics
######################################################################

# load the data
library(qtl)
library(qtlbook)
data(ovar)

# summary of the data
summary(ovar)

# summary plot of the data
plot(ovar)

# scatterplot of on1 and on2
plot(jitter(on2) ~ jitter(on1), data=ovar$pheno,
     xlab="on1", ylab="on2", cex=0.6,
     xlim=c(6.66, 17.53), ylim=c(6.66, 17.53))

# correlation between on1 and on2
cor(pull.pheno(ovar,"on1"), pull.pheno(ovar,"on2"), 
    use="complete")

# check that onm = (on1+on2)/2
max(abs(pull.pheno(ovar,"onm") - 
        (pull.pheno(ovar,"on1")+pull.pheno(ovar,"on2"))/2),
    na.rm=TRUE)

# check whether onm is missing when on1 or on2 are
table(apply(is.na(pull.pheno(ovar,1:3)), 1, paste,
            collapse=":"))

# individuals with onm but not on1
ovar$pheno[!is.na(ovar$pheno$onm) & is.na(ovar$pheno$on1),]

# Make onm phenotype missing when on1 is missing
ovar$pheno$onm[is.na(ovar$pheno$on1)] <- NA

# box plot of ovariole counts in the two crosses
boxplot(onm ~ cross, data=ovar$pheno, horizontal=TRUE,
        xlab="Ovariole count", ylab="Cross")

# t-test for difference in average phenotype between the two crosses
t.test(onm ~ cross, data=ovar$pheno)

# pull out first cross
ovar1 <- subset(ovar, ind=(pull.pheno(ovar, "cross")==1))

# plot phenotype against number of typed genotypes
plot(jitter(ntyped(ovar1)),jitter(pull.pheno(ovar1,"onm")),
     xlab="No. genotypes", ylab="Ovariole count")

# plot pairwise recombination fractions
ovar1 <- est.rf(ovar1)
plot.rf(ovar1)

# reestimate map and plot it against that in data
newmap <- est.map(ovar1, error.prob=0.001, verbose=FALSE)
plot.map(ovar1, newmap)

# estimate map with Kosambi map function
newmap.k <- est.map(ovar1, err=0.001, map.function="kosambi")

# summary of the two maps
summary(newmap)
summary(newmap.k)

# two-locus genotypes for markers st and e in 2nd cross
ovar2 <- subset(ovar, ind=(pull.pheno(ovar,"cross")==2))
geno.crosstab(ovar2, "st", "e")

# two-locus genotypes for markers st and e in 1st cross
geno.crosstab(ovar1, "st", "e")

# plot chr 3 genotypes for a random set of individuals from the 2nd cross
toplot <- sort(sample(nind(ovar2), 15))
plot.geno(ovar2, chr=3, ind=toplot)

######################################################################
# 10.2 Initial cross
######################################################################

# Omit individuals with missing genotype data
ovar1 <- subset(ovar1, ind=!is.na(pull.pheno(ovar1,"onm")))

# Multiple imputations
ovar1 <- sim.geno(ovar1, n.draws=512, step=2, err=0.001, 
                  map.function="kosambi")

# genome scan by multiple imputation
out1 <- scanone(ovar1, method="imp")

# summary of genome scan
summary(out1)

# plot of genome scan
plot(out1, ylab="LOD score")

# effect plot for chr 3
effectplot(ovar1, mname1="cpo")

# fitqtl to estimate QTL effect
qtl <- makeqtl(ovar1, chr=3, pos=82.4)
summary(fitqtl(ovar1, qtl=qtl, dropone=FALSE, get.ests=TRUE))

# scan adjusting for chr 3 locus
out1.c3 <- addqtl(ovar1, qtl=qtl)

# summary of genome scan, adjusting for chr 3 locus
summary(out1.c3)

# plot of genome scan, adjusting for chr 3 locus
plot(out1.c3, ylab="LOD score")

# 2d scan on chromosome 3
qtl2 <- makeqtl(ovar1, 2, 90)
out1.ap <- addpair(ovar1, qtl=qtl2, chr=3, verbose=FALSE)

# summary of the 2d scan
summary(out1.ap)

# refine QTL locations in the 3-QTL model
qtl <- makeqtl(ovar1, c(2, 3, 3), c(90, 82.8, 121))
rqtl <- refineqtl(ovar1, qtl=qtl, verbose=FALSE)

# print the results of refineqtl
rqtl

# plot LOD profiles for the 3-QTL model
plotLodProfile(rqtl, col=c("black","red","blue"), 
               ylab="Profile LOD score")

# approximate confidence intervals for the location of the proximal c3 QTL
lodint(rqtl, qtl.index=2)
bayesint(rqtl, qtl.index=2)

# permutation test with 2-dim, 2-QTL genome scan
strat <- (nmissing(ovar1) < 15)
operm <- scantwo(ovar1, method="imp", n.perm=1000, 
                 perm.strata=strat)

# Significance thresholds
summary(operm)

# Calculate penalties 
print(pen <- calc.penalties(operm))

# stepwise selection starting at a 3-qtl model
stepout1 <- stepwiseqtl(ovar1, penalties=pen, qtl=rqtl, 
                        max.qtl=8, verbose=FALSE)

# Print stepwiseqtl results
stepout1

# Plots to study the inferred pair of QTL on chr 2
mar <- find.marker(ovar1, 2, c(92.4, 94))
par(mfrow=c(1,2))
effectplot(ovar1, mname1="2@92.4", mname2="2@94")
plot.pxg(ovar1, marker=mar)

# drop a QTL, refine locations, do dropone analysis
qtl <- dropfromqtl(stepout1, 2)
qtl <- refineqtl(ovar1, qtl=qtl, verbose=FALSE)

# Print the result
qtl

######################################################################
# 10.3 Combined data
######################################################################

# clean up memory
object.size(ovar1)

object.size(ovar1)/1024^2

ovar1 <- clean(ovar1)
object.size(ovar1)

object.size(ovar1)/1024

# remove individuals with missing phenotype
ovar <- subset(ovar, ind=!is.na(pull.pheno(ovar,"onm")))

# Multiple imputations
ovar <- sim.geno(ovar, n.draws=512, step=2, err=0.001, 
                 map.function="kosambi")

# genome scan with combined data
cross <- pull.pheno(ovar, "cross")
outc <- scanone(ovar, method="imp", addcovar=cross)

# genome scan with just the second cross
out2 <- scanone(subset(ovar, ind=(cross==2)), method="imp")

# plot genome scan for combined data and the two individual crosses
plot(outc, out1, out2, ylab="LOD score")

# location of maximum LOD in combined data
max(outc)

cross <- data.frame(cross=cross)
qtlc <- makeqtl(ovar, chr=3, pos=72.9)
outc.c3 <- addqtl(ovar, qtl=qtlc, covar=cross)
qtl2 <- makeqtl(subset(ovar, ind=(cross==2)), chr=3, pos=72.9)
out2.c3 <- addqtl(subset(ovar, ind=(cross==2)), qtl=qtl2)

# plot genome scan for combined data and the two individual crosses
plot(outc.c3, out1.c3, out2.c3, ylab="LOD score")

# create vector indicating the 3 strata
strat <- pull.pheno(ovar, "cross")
strat[strat==1 & nmissing(ovar)<15] <- 3
table(strat)

# the permutation test
opermc <- scantwo(ovar, method="imp", addcovar=cross,
                  n.perm=n.perm, perm.strata=strat)

# significance thresholds
summary(opermc)
print(penc <- calc.penalties(opermc))

# 2d scan on chr 3, controlling for chr 2 locus
qtl.c2 <- makeqtl(ovar, 2, 115)
out.ap.c3 <- addpair(ovar, qtl=qtl.c2, chr=3, covar=cross, 
                   verbose=FALSE)

# summary of the results
summary(out.ap.c3)

# plot of the results
plot(out.ap.c3, lower="cond-int")

# 2d scan on chr 2, controlling for the chr 3 locus
out.ap.c2 <- addpair(ovar, qtl=qtlc, chr=2, covar=cross, 
                   verbose=FALSE)

# summary of the results
summary(out.ap.c2, perms=opermc, pval=TRUE)

# plot of the results
plot(out.ap.c2, lower="cond-int", upper="cond-add")

# refine qtl locations
qtl <- makeqtl(ovar, c(2, 2, 3, 3), c(4, 114, 62.8, 74.8))
rqtl <- refineqtl(ovar, qtl=qtl, covar=cross, verbose=FALSE,
                  formula=y~cross+Q1+Q2+Q3+Q4+Q3:Q4)

# print result
rqtl

# drop-one analysis with fitqtl
summary(fitqtl(ovar, qtl=rqtl, covar=cross, 
               formula=y~cross+Q1+Q2+Q3+Q4+Q3:Q4), 
        pvalues=FALSE)

# plot lod profiles
plotLodProfile(rqtl, col=c("red","blue","red","blue"),
               ylab="Profile LOD score")

addint(ovar, qtl=rqtl, covar=cross, 
       formula=y~cross+Q1+Q2+Q3+Q4+Q3:Q4,
       pvalues=FALSE)

# scan for one more QTL
onemore <- addqtl(ovar, qtl=rqtl, covar=cross, 
                  formula=y~cross+Q1+Q2+Q3+Q4+Q3:Q4)

# maximum LOD on each chromosome
summary(onemore)

qtl2 <- addtoqtl(ovar, rqtl, 2, 89.1)
qtl2 <- reorderqtl(qtl2)
rqtl2 <- refineqtl(ovar, qtl=qtl2, covar=cross, verbose=FALSE,
                   formula=y~cross+Q1+Q2+Q3+Q4+Q5+Q4:Q5)

# print result
rqtl2

summary(fitqtl(ovar, qtl=rqtl2, covar=cross, 
               formula=y~cross+Q1+Q2+Q3+Q4+Q5+Q4:Q5),
        pvalues=FALSE)

# stepwise QTL analysis 
stepout1 <- stepwiseqtl(ovar, covar=cross, pen=pen,
                        max.qtl=8, verbose=FALSE)

# print results
stepout1

# stepwise QTL analysis with exclusively heavy penalties
stepout2 <- stepwiseqtl(ovar, covar=cross, pen=pen[1:2],
                        max.qtl=8, verbose=FALSE)

# print results
stepout2

######################################################################
# 10.4 Discussion
######################################################################

# end of chap10.R
