##############################################################
# R code for "R/qtl tutorial"
#
# Karl W Broman, broman@wisc.edu
# University of Wisconsin Madison
#
# https://rqtl.org
#
# 30 November 2012
##############################################################

############################################################
# Preliminaries
############################################################
library(qtl)

ls()

help(read.cross)
?read.cross

# url.show("https://rqtl.org/rqtltour3.R")

############################################################
# Data import
############################################################
sug <- read.cross("csv", "https://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

############################################################
# Summaries
############################################################
summary(sug)

nind(sug)
nchr(sug)
totmar(sug)
nmar(sug)
nphe(sug)

plot(sug)

plotMissing(sug)
plotMap(sug)
plotPheno(sug, pheno.col=1)
plotPheno(sug, pheno.col=2)
plotPheno(sug, pheno.col=3)
plotPheno(sug, pheno.col=4)
plotPheno(sug, pheno.col=5)
plotPheno(sug, pheno.col=6)

plotPheno(sug, pheno.col="bp")
plotPheno(sug, pheno.col="bw")

############################################################
# Single-QTL analysis
############################################################
sug <- calc.genoprob(sug, step=1)

out.em <- scanone(sug)

summary(out.em)

summary(out.em, threshold=3)

plot(out.em)

out.hk <- scanone(sug, method="hk")

plot(out.em, out.hk, col=c("blue", "red"))

plot(out.em, col="blue")
plot(out.hk, col="red", add=TRUE)

plot(out.em, out.hk, col=c("blue", "red"), chr=c(7,15))

plot(out.em, col="blue", chr=c(7,15))
plot(out.hk, col="red", chr=c(7,15), add=TRUE)

plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")

############################################################
# Permutation tests
############################################################
load(url("https://rqtl.org/various.RData"))

operm <- scanone(sug, method="hk", n.perm=1000)

plot(operm)

summary(operm)
summary(operm, alpha=c(0.05, 0.2))

summary(out.hk, perms=operm, alpha=0.2, pvalues=TRUE)

############################################################
# Interval estimates of QTL location
############################################################
lodint(out.hk, chr=7)
bayesint(out.hk, chr=7)

lodint(out.hk, chr=7, expandtomarkers=TRUE)
bayesint(out.hk, chr=7, expandtomarkers=TRUE)

lodint(out.hk, chr=7, drop=2)
bayesint(out.hk, chr=7, prob=0.99)

lodint(out.hk, chr=15)
bayesint(out.hk, chr=15)

############################################################
# QTL effects
############################################################
max(out.hk)
mar <- find.marker(sug, chr=7, pos=47.7)
plotPXG(sug, marker=mar)

sug <- sim.geno(sug, n.draws=64, step=1)
effectplot(sug, mname1=mar)

effectplot(sug, mname1="7@47.7")

max(out.hk, chr=15)
mar2 <- find.marker(sug, chr=15, pos=12)
plotPXG(sug, marker=mar2)
effectplot(sug, mname1="15@12")

plotPXG(sug, marker=c(mar, mar2))
plotPXG(sug, marker=c(mar2, mar))

effectplot(sug, mname1="7@47.7", mname2="15@12")
effectplot(sug, mname2="7@47.7", mname1="15@12")

############################################################
# Other phenotypes
############################################################
out.hr <- scanone(sug, pheno.col=2, method="hk")

out.bw <- scanone(sug, pheno.col="bw", method="hk")

out.logbw <- scanone(sug, pheno.col=log(sug$pheno$bw), method="hk")

out.all <- scanone(sug, pheno.col=1:4, method="hk")

summary(out.all, threshold=3)

summary(out.all, threshold=3, lodcolumn=4)

summary(out.all, threshold=3, format="allpeaks")

summary(out.all, threshold=3, format="allpheno")

summary(out.all, threshold=3, format="tabByCol")
summary(out.all, threshold=3, format="tabByChr")

############################################################
# Data diagnostics
############################################################
bad <- read.cross("csv", "https://rqtl.org", "bad.csv")

summary(bad)

par(mfrow=c(3,1))
for(i in 1:3)
  plotPheno(bad, pheno.col=i)

par(mfrow=c(3,1))
for(i in 1:3)
  plot(bad$pheno[,i], ylab=names(bad$pheno)[i])

par(mfrow=c(1,1))
pairs(bad$pheno)

par(mfrow=c(2,1))
plot(ntyped(bad, "ind"), main="No. genotypes, by ind'l",
            ylab="No. genotypes")
plot(ntyped(bad, "mar"), main="No. genotypes, by marker",
            ylab="No. genotypes")

bad <- subset(bad, ind=(ntyped(bad, "ind") > 100))
nt.mar <- ntyped(bad, "mar")
bad <- drop.markers(bad, names(nt.mar)[nt.mar < 180])

cg <- comparegeno(bad)
lowertri <- cg[lower.tri(cg)]
par(mfrow=c(1,1))
hist(lowertri, breaks=seq(0, 1, by=0.01))
rug(lowertri)

sort(lowertri, decreasing=TRUE)[1:10]

wh <- which(!is.na(cg) & cg > 0.9, arr.ind=TRUE)
bad <- subset(bad, ind = -wh[,2])

gt <- geno.table(bad)
gt[gt$P.value < 0.05/nrow(gt),]

gt2 <- geno.table(bad, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt2)
plot(gt2, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="green")
legend("bottomleft", colnames(gt2)[5:7], lwd=2,
       col=c("black","blue","red"))

geno.table(bad, chr=3)

geno.table(bad, chr=6)

gt3 <- geno.table(bad, chr=-3)
badmar <- rownames(gt3)[gt3$P.value < 0.05/nrow(gt3)]
bad <- drop.markers(bad, badmar)

bad <- est.rf(bad)
par(mfrow=c(1,1))
plotRF(bad)

plotRF(bad, chr=c(7,10))

markernames(bad, chr=7)
markernames(bad, chr=7)[4]

newmap <- est.map(bad)
par(mfrow=c(1,1))
plotMap(bad, newmap)

plotMap(bad, newmap, chr=7, show.marker.names=TRUE)

out <- tryallpositions(bad, "D7M4", chr=10)
summary(out)

bad <- movemarker(bad, "D7M4", 10, 9.245)

newmap <- est.map(bad)
plotMap(bad, newmap)

plotRF(bad, chr=18)

rip <- ripple(bad, chr=18, window=7)
summary(rip)

compareorder(bad, chr=18, rip[2,])

mar <- markernames(bad, chr=18)[3]
out <- tryallpositions(bad, mar, 18)
summary(out)

bad <- movemarker(bad, mar, 18, 53.7)
newmap <- est.map(bad)
plotMap(bad, newmap)

xo <- countXO(bad)
hist(xo, breaks=50)
rug(xo)

bad <- subset(bad, ind=(xo < 60))

bad <- calc.errorlod(bad, err=0.01)
top.errorlod(bad)

############################################################
# Two-dimensional, two-QTL scans
############################################################
sug <- read.cross("csv", "https://rqtl.org", "sug.csv",
                  genotypes=c("CC", "CB", "BB"), alleles=c("C", "B"))

load(url("https://rqtl.org/various.RData"))

sug <- calc.genoprob(sug, step=2)

# out2 <- scantwo(sug, method="hk")

plot(out2)

plot(out2, lower="fv1")

plot(out2, lower="fv1", upper="av1")

# operm2 <- scantwo(sug, method="hk", n.perm=5)

summary(out2, perms=operm2, alpha=0.2, pvalues=TRUE)

############################################################
# Multiple-QTL analyses
############################################################
sug <- calc.genoprob(sug, step=1)

qtl <- makeqtl(sug, chr=c(7,15), pos=c(47.7, 12), what="prob")

out.fq <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq)

summary(fitqtl(sug, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))

out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1*Q2)
out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1+Q2+Q1:Q2)
summary(out.fqi)

addint(sug, qtl=qtl, method="hk")

rqtl <- refineqtl(sug, qtl=qtl, method="hk")
rqtl

summary(out.fqr <- fitqtl(sug, qtl=rqtl, method="hk"))

plotLodProfile(rqtl)

out.hk <- scanone(sug, method="hk")
plot(out.hk, chr=c(7,15), col="red", add=TRUE)

out.aq <- addqtl(sug, qtl=rqtl, method="hk")

plot(out.aq)

print(pen <- calc.penalties(operm2))

out.sq <- stepwiseqtl(sug, max.qtl=5, penalties=pen, method="hk", verbose=2)
out.sq

# end of rqtltour3.R
