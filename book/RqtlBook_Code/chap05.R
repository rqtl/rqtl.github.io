######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 5: Non-normal phenotypes
#######################################################################

######################################################################
# 5.1 Nonparametric interval mapping
######################################################################

# Get access to the listeria data and run calc.genoprob
library(qtl)
data(listeria)
listeria <- calc.genoprob(listeria, step=1, error.prob=0.001)

# Nonparametric interval mapping
out.np <- scanone(listeria, model="np")

# Plot results of scanone with method="np"
plot(out.np, ylab="LOD score", alternate.chrid=TRUE)

# Permutation test for nonparametric interval mapping
operm.np <- scanone(listeria, model="np", n.perm=1000, 
                    perm.Xsp=TRUE)

# LOD threshold for nonparametric interval mapping
summary(operm.np, 0.05)

# Summary of the significant loci
summary(out.np, perms=operm.np, alpha=0.05, pvalues=TRUE)

######################################################################
# 5.2 Binary traits
######################################################################

# Create a binary phenotype for the listeria data
binphe <- as.numeric(pull.pheno(listeria, 1) > 250)
listeria$pheno <- cbind(listeria$pheno, binary=binphe)

# Binary trait mapping
out.bin <- scanone(listeria, pheno.col="binary", 
                   model="binary")

# Plot results of scanone with method="binary" and method="np"
plot(out.np, out.bin, col=c("blue", "red"), ylab="LOD score",
     alternate.chrid=TRUE)

# Permutation test for binary trait mapping
operm.bin <- scanone(listeria, pheno.col="binary", 
                     model="binary", n.perm=1000, 
                     perm.Xsp=TRUE)

# LOD threshold for binary trait mapping
summary(operm.bin, alpha=0.05)

# Summary of the significant loci
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE)

######################################################################
# 5.3 Two-part model
######################################################################

# Attach log(survival) to the phenotype data
y <- log(pull.pheno(listeria, 1))
listeria$pheno <- cbind(listeria$pheno, logsurv=y)

# Interval mapping with the two-part model
out.2p <- scanone(listeria, model="2part", upper=TRUE, 
                  pheno.col="logsurv")

# Plot the two-part model results
plot(out.2p, lodcolumn=1:3, ylab="LOD score", 
     alternate.chrid=TRUE)

# Permutation test for the two-part model
operm.2p <- scanone(listeria, model="2part", upper=TRUE,
                    pheno.col="logsurv", n.perm=1000, 
                    perm.Xsp=TRUE)

# LOD thresholds by 2-part model
summary(operm.2p, alpha=0.05)

# Summary of the significant loci
summary(out.2p, perms=operm.2p, alpha=0.05, pvalues=TRUE)

# Alternate summary
summary(out.2p, perms=operm.2p, alpha=0.05, pvalues=TRUE,
        format="allpeaks")

######################################################################
# 5.4 Other extensions
######################################################################

# Need the survival package and the scanone_cph.R file
library(survival)
source("scanone_cph.R")

# Create a phenotype that is a censored survival time
y <- pull.pheno(listeria, 1)
y <- Surv(y, y<250)
listeria$pheno <- cbind(listeria$pheno, surv=y)

# Haley--Knott regr with Cox prop'n hazards model
out.cph <- scanone.cph(listeria, pheno.col="surv")

# Plot the Cox proportional hazards model results
plot(out.cph, ylab="LOD score", alternate.chrid=TRUE)

# permutation test
n.perm <- 1000
operm.cph <- cbind(lod=1:n.perm)
chr <- names(listeria$geno)
temp <- subset(listeria, chr=(chr != "X"))
n.ind <- nind(listeria)
for(i in 1:n.perm) {
  temp$pheno <- temp$pheno[sample(n.ind),]
  out <- scanone.cph(temp, pheno.col="surv")
  operm.cph[i] <- max(out[,3], na.rm=TRUE)
}
class(operm.cph) <- "scanoneperm"

# LOD threshold for Cox proportional hazards model
summary(operm.cph, 0.05)

# Summary of the significant loci
summary(out.cph, perms=operm.cph, alpha=0.05, pvalues=TRUE)

# end of chap05.R
