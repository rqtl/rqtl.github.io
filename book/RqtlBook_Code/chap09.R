######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 9: Fit and exploration of multiple-QTL models
#######################################################################

######################################################################
# 9.1 Model selection
######################################################################

###################################
# 9.1.1 Class of models
###################################

###################################
# 9.1.2 Model fit
###################################

###################################
# 9.1.3 Model search
###################################

###################################
# 9.1.4 Model comparison
###################################

###################################
# 9.1.5 Further discussion
###################################

######################################################################
# 9.2 Bayesian QTL mapping
######################################################################

######################################################################
# 9.3 Multiple QTL mapping in R/qtl
######################################################################

###################################
# 9.3.1 makeqtl and fitqtl
###################################

# load data
library(qtl)
data(hyper)

# Perform multiple imputations
hyper <- sim.geno(hyper, step=2, n.draws=128, err=0.001)

# Create QTL object
qtl <- makeqtl(hyper, chr=c(1, 4, 6, 15), 
               pos=c(68.3, 30, 60, 18))

# print QTL object
qtl

# Plot QTL object
plot(qtl)

# Fit multiple QTL model
out.fq <- fitqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4)
summary(out.fq)

# Get estimated QTL effects for multiple QTL model
out.fq2 <- fitqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4,
                  dropone=FALSE, get.ests=TRUE)
summary(out.fq2)

###################################
# 9.3.2 refineqtl
###################################

# Refine estimates of QTL locations 
rqtl <- refineqtl(hyper, qtl=qtl, formula=y~Q1+Q2+Q3*Q4, 
                  verbose=FALSE)

# Print refined QTL object
rqtl

# Fit multiple QTL object to see improvement in model fit with new locations
out.fq3 <- fitqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4, 
                  dropone=FALSE)
summary(out.fq3)

# plot LOD profile
plotLodProfile(rqtl, ylab="Profile LOD score")

# approximate confidence intervals for the location of the QTL on chr 4
lodint(rqtl, qtl.index=2)
bayesint(rqtl, qtl.index=2)

###################################
# 9.3.3 addint
###################################

# Study the addition of each possible pairwise interaction
addint(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4, pvalues=FALSE)

# Study addition of a pairwise interaction if the 6x15 interaction is not included
addint(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3+Q4, pvalues=FALSE)

###################################
# 9.3.4 addqtl
###################################

# Scan for an additional additive QTL
out.aq <- addqtl(hyper, qtl=rqtl, formula=y~Q1+Q2+Q3*Q4)

# Maximum LOD peak
max(out.aq)

# Plot the LOD curves
plot(out.aq, ylab="LOD score")

# Scan for an additional QTL, interacting with the chr 15 locus
out.aqi <- addqtl(hyper, qtl=rqtl, 
                  formula=y~Q1+Q2+Q3*Q4+Q4*Q5)

# Plot the LOD curves
plot(out.aqi, ylab="LOD score")

# Plot the epistasis LOD curves
plot(out.aqi - out.aq, ylab="LOD score")

###################################
# 9.3.5 addpair
###################################

# Scan for a pair of QTL on chr 1
out.ap <- addpair(hyper, qtl=rqtl, chr=1, formula=y~Q2+Q3*Q4, 
                  verbose=FALSE)

# Look at evidence for a second QTL
summary(out.ap)

# Plot the 2d scan results
plot(out.ap, lower="cond-int", upper="cond-add")

# Scan for an interacting pair, with one locus interacting with the chr 6 QTL
out.ap2 <- addpair(hyper, qtl=rqtl, chr=c(7,15), verbose=TRUE,
                   formula=y~Q1+Q2+Q3+Q5*Q6+Q3:Q5)

# Summarize the 2d scan results
summary(out.ap2)

# Plot the 2d scan results
plot(out.ap2)

###################################
# 9.3.6 Manipulating qtl objects
###################################

# Add an additional QTL to a QTL object
print( rqtl2 <- addtoqtl(hyper, rqtl, 1, 43.3) )

# Move one QTL to a different location
print( rqtl3 <- replaceqtl(hyper, rqtl2, 1, 1, 77.3) )

# Reorder the QTL in a QTL object
print( rqtl4 <- reorderqtl(rqtl3, c(4,3,2,1,5)) )

# Reorder the QTL in a QTL object according to genomic positions
print( rqtl5 <- reorderqtl(rqtl4) )

# Drop a QTL from the object
print( rqtl6 <- dropfromqtl(rqtl5, 1) )

###################################
# 9.3.7 stepwiseqtl
###################################

# Stratified permutation test for 2d, 2-QTL scan
strat <- (nmissing(hyper) > 50)
operm2 <- scantwo(hyper, method="imp", n.perm=1000, 
                  perm.strat=strat)

# Estimated significance thresholds
summary(operm2)

# Calculate penalties for use with stepwiseqtl
print(pen <- calc.penalties(operm2))

# Calculate penalties for both 5% and 20% significance levels
calc.penalties(operm2, alpha=c(0.05, 0.20))

# stepwise model selection 
outsw1 <- stepwiseqtl(hyper, max.qtl=8, penalties=pen, 
                      verbose=FALSE)

# print results
outsw1

# stepwise model selection, saving LOD profiles and trace through model space
outsw2 <- stepwiseqtl(hyper, max.qtl=8, penalties=pen, 
                      verbose=FALSE, keeplodprofile=TRUE,
                      keeptrace=TRUE)

# plot LOD profiles
plotLodProfile(outsw2)

# print attribute names
names(attributes(outsw2))

# pull out the model "trace" and print the first model seen
thetrace <- attr(outsw2, "trace")
thetrace[[1]]

# plot models visited in stepwiseqtl
par(mfrow=c(6,3))
for(i in seq(along=thetrace))
  plotModel(thetrace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =", 
              round(attr(thetrace[[i]], "pLOD"), 2)))

# stepwise analysis with all heavy penalties
outsw3 <- stepwiseqtl(hyper, max.qtl=8, penalties=pen[1:2],
                      verbose=FALSE) 

# print results
outsw3

# stepwise analysis with more liberal penalties
liberalpen <- calc.penalties(operm2, alpha=0.2)
outsw4 <- stepwiseqtl(hyper, max.qtl=8, penalties=liberalpen,
                      verbose=FALSE) 

# print results
outsw4

# stepwise analysis allowing only additive models
outsw5 <- stepwiseqtl(hyper, max.qtl=8, penalties=pen,
                      additive.only=TRUE, verbose=FALSE)

# print results
outsw5

# stepwise analysis starting at the 4-QTL model with the 6x15 interaction
qtl <- makeqtl(hyper, chr=c(1, 4, 6, 15), 
               pos=c(68.3, 30, 60, 18))
outsw6 <- stepwiseqtl(hyper, max.qtl=8, penalties=pen,
                      qtl=qtl, formula=y~Q1+Q2+Q3*Q4, 
                      verbose=FALSE)

# print results
outsw6

# end of chap09.R
