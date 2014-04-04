######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 8: Two-dimensional, two-QTL scans
#######################################################################

######################################################################
# 8.1 The normal model
######################################################################

# load hyper and run calc.genoprob
library(qtl)
data(hyper)
hyper <- calc.genoprob(hyper, step=2.5, err=0.001)

# two-dimensional scan
out2 <- scantwo(hyper, verbose=FALSE)

# Plot the results
plot(out2, chr=c(1,4,6,7,15))

# Plot the results with LOD_fv1 in the lower right triangle
plot(out2, chr=c(1,4,6,7,15), lower="cond-int")

# Plot the results with LOD_fv1 in the lower right triangle
plot(out2, chr=c(1,4,6,7,15), upper="add", lower="cond-add")

# Plot the results for chromosome 1
plot(out2, chr=1, lower="cond-int", upper="cond-add")

# Summary of 2d scan
summary(out2, thresholds=c(6.0, 4.7, 4.4, 4.7, 2.6))

# Summary of 2d scan, ignored M_i by taking T_i = Inf
summary(out2, thresholds=c(6.0, 4.7, Inf, 4.7, 2.6))

# Stratified permutation test for the 2d scan
strat <- (nmissing(hyper) > 50)
operm2 <- scantwo(hyper, n.perm=1000, perm.strat=strat)

# Thresholds for the 2d scan results 
summary(operm2)

# Permutation test in two batches
set.seed(85842518)
operm2a <- scantwo(hyper, n.perm=500, perm.strat=strat)
save(operm2a, file="perm2a.RData")

# Second batch
set.seed(85842519)
operm2b <- scantwo(hyper, n.perm=500, perm.strat=strat)
save(operm2b, file="perm2b.RData")

# Combine the batches
load("perm2a.RData")
load("perm2b.RData")
operm2 <- c(operm2a, operm2b)

# summary of 2d scan results, using the permutation results and including p-values
summary(out2, perms=operm2, alphas=c(0.2, 0.2, 0, 0.2, 0.2),
        pvalues=TRUE)

# Plot the results for chromosome 3, zeroing out bits along the diagonal
plot(clean(out2), chr=3, lower="cond-add")

# summary of the "cleaned" 2d scan results
summary(clean(out2), perms=operm2, 
        alphas=c(0.2, 0.2, 0, 0.2, 0.2))

# Impute genotypes on chr 3
hyperc3 <- sim.geno(subset(hyper, chr=3), step=2.5,
                    error.prob=0.001, n.draws=256)

# find markers near the inferred QTL
mar <- find.marker(hyperc3, "3", c(37.2, 44.7))

# Plot of effects for tightly linked loci on chr 3
par(mfrow=c(1,2))
plot.pxg(hyperc3, marker=mar)
effectplot(hyperc3, mname1="3@37.2", mname2="3@44.7",
           ylim=range(pull.pheno(hyperc3,1)))

# Plot of effects for tightly linked loci on chr 3
hypersub <- sim.geno(subset(hyper, chr=c(1,4,6,15)), step=2.5,
                     error.prob=0.001, n.draws=256)
par(mfrow=c(1,2))
effectplot(hypersub, mname1="1@68.3", mname2="4@30",
           ylim=c(95, 110))
effectplot(hypersub, mname1="6@60", mname2="15@18",
           ylim=c(95, 110))

######################################################################
# 8.2 Binary traits
######################################################################

# load the nf1 data
library(qtlbook)
data(nf1)
nf1 <- drop.nullmarkers(nf1)

# Calculate genotype probabilities
nf1 <- calc.genoprob(nf1, step=2.5, error.prob=0.001)
frommom <- pull.pheno(nf1,"frommom")

# Perform the two-QTL analyses
out.frommom <- scantwo(subset(nf1, ind=(frommom==1)), 
                       model="binary")
out.fromdad <- scantwo(subset(nf1, ind=(frommom==0)), 
                       model="binary")

# Perform permutation test for the two-QTL analyses
operm.frommom <- scantwo(subset(nf1, ind=frommom==1), 
                         model="binary", n.perm=1000)
operm.fromdad <- scantwo(subset(nf1, ind=frommom==0), 
                         model="binary", n.perm=1000)

# Thresholds from permutations
summary(operm.frommom)
summary(operm.fromdad)

# Summary of 2d scan results for mutation from mom
summary(out.frommom, perms=operm.frommom, pvalues=TRUE,
        alphas=c(0.2, 0.2, 0, 0.2, 0.2))

# Plot the results for 2d scan, mutation from mom 
plot(out.frommom, chr=c(7,15,17), lower="cond-int")

# Summary of 2d scan results for mutation from dad
summary(out.fromdad, perms=operm.fromdad, pvalues=TRUE,
        alphas=c(0.2, 0.2, 0, 0.2, 0.2))

# Plot the results for 2d scan, mutation from dad 
plot(out.fromdad, chr=c(9, 19), lower="cond-add")

# Preparation for effectplot
nf1.fm <- sim.geno(subset(nf1, chr=c(7,17), ind=(frommom==1)),
                   step=2.5, error.prob=0.001, n.draws=256)
nf1.fd <- sim.geno(subset(nf1, chr=c(9,19), ind=(frommom==0)),
                   step=2.5, error.prob=0.001, n.draws=256)

# Plot estimated proportions of affected individuals vs QTL genotypes
par(mfrow=c(1,2))
effectplot(nf1.fm, mname1="7@45", mname2="17@3",
           ylim=c(0,1))
effectplot(nf1.fd, mname1="9@55.5", mname2="19@0",
           ylim=c(0,1))

######################################################################
# 8.3 The X chromosome
######################################################################

# Calculate genotype probabilities
data(gutlength)
gutlength <- calc.genoprob(gutlength, step=5,
                           error.prob=0.001)

# 2d scan for gutlength data
out.gl <- scantwo(gutlength)

# Plot the results for 2d scan
plot(out.gl, lower="cond-int", alternate.chrid=TRUE)

# Significant peaks in 2d scan
summary(out.gl, thresholds=c(9.1, 7.1, 6.3, 6.3, 3.3))

######################################################################
# 8.4 Covariates
######################################################################

# form covariate matrix
cross <- as.numeric(pull.pheno(gutlength,"cross"))
frommom <- as.numeric(cross < 3)
forw <- as.numeric(cross == 1 | cross == 3)
sex <- as.numeric(pull.pheno(gutlength,"sex") == "M")
crossX <- cbind(frommom, forw, frommom*forw)
x <- cbind(sex, crossX)

# 2d scan with additive covariates
out.gl.a <- scantwo(gutlength, addcovar=x)

# Plot the differences in LOD scores, with and without the covariates
plot(out.gl.a - out.gl, allow.neg=TRUE, alternate.chrid=TRUE)

# maximum LOD scores with and without covariates
max(out.gl)
max(out.gl.a)

# end of chap08.R
