######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 6: Experimental design and power
#######################################################################

######################################################################
# 6.1 Phenotypes and covariates
######################################################################

######################################################################
# 6.2 Strains and strain surveys
######################################################################

######################################################################
# 6.3 Theory
######################################################################

###################################
# 6.3.1 Variance attributable to a locus
###################################

###################################
# 6.3.2 Residual error variance
###################################

###################################
# 6.3.3 Information content
###################################

######################################################################
# 6.4 Examples with R/qtlDesign
######################################################################

###################################
# 6.4.1 Functions
###################################

###################################
# 6.4.2 Choosing a cross
###################################

# Load R/qtlDesign
library(qtlDesign)

# LOD threshold for a mouse backcross
thresh(G=1440, cross="bc", p=0.05)

# LOD threshold for a mouse intercross
thresh(G=1440, cross="f2", p=0.05)

# Detectable effects for a sample size of 100
detectable(cross="bc", n=100, sigma2=1, thresh=3.2)
detectable(cross="f2", n=100, sigma2=1, thresh=4.2)

# Sample size in intercross
samplesize(cross="f2", effect=c(5,0), env.var=64, gen.var=25,
           thresh=4.2)

# Sample size in intercross for larger genetic variance
samplesize(cross="f2", effect=c(5,0), env.var=64, gen.var=50,
           thresh=4.2)
samplesize(cross="f2", effect=c(5,0), env.var=64, 
           gen.var=100, thresh=4.2)

# Sample sizes in a backcross
samplesize(cross="bc", effect=5, env.var=64, gen.var=25,
           thresh=3.2)
samplesize(cross="bc", effect=5, env.var=64, gen.var=50,
           thresh=3.2)
samplesize(cross="bc", effect=5, env.var=64, gen.var=100,
           thresh=3.2)

# LOD threshold in RIL by sibling mating
thresh(G=1440*4, cross="bc", p=0.05)

# Sample size in RIL
samplesize(cross="ri", effect=5, env.var=64, gen.var=25,
           thresh=3.8)

# RIL sample size with biological replicates
samplesize(cross="ri", effect=5, env.var=64, gen.var=25,
           thresh=3.8, bio.rep=2)
samplesize(cross="ri", effect=5, env.var=64, gen.var=25,
           thresh=3.8, bio.rep=4)
samplesize(cross="ri", effect=5, env.var=64, gen.var=25,
           thresh=3.8, bio.rep=16)
samplesize(cross="ri", effect=5, env.var=64, gen.var=25,
           thresh=3.8, bio.rep=100)

###################################
# 6.4.3 Genotyping strategies
###################################

# Optimal marker spacing for intercross when all individuals are typed
optspacing(cost=1/300, G=1440, sel.frac=1, cross="f2")

# Optimal marker spacing and selection fraction for an intercross
optspacing(cost=10/3000, G=1440, sel.frac=NULL, cross="f2")

###################################
# 6.4.4 Phenotyping strategies
###################################

# Plot illustrating phenotyping strategies
library(qtl)
mp <- sim.map(len=rep(100,5), n.mar=11, include.x=FALSE,
              eq.spacing=TRUE)
cr <- sim.cross(mp, model=c(1,50,1,0), n.ind=200, type="f2")
idx40 <- mma(pull.geno(cr, chr=1), p=40)
cr <- calc.genoprob(cr, step=2)
out1 <- scanone(cr)
out2 <- scanone(subset(cr, ind=idx40$cList))
out3 <- scanone(subset(cr, ind=sample(1:200, 40)))
plot(out1, out2, out3, ylab="LOD score")

###################################
# 6.4.5 Fine mapping
###################################

# Median length of confidence intervals in backcross and intercross
ci.length(cross="bc", n=250, effect=5, p=0.95, 
          gen.var=25, env.var=64)
ci.length(cross="f2", n=250, effect=c(5,0), p=0.95, 
          gen.var=25, env.var=64)

######################################################################
# 6.5 Other experimental populations
######################################################################

######################################################################
# 6.6 Estimating power and precision by simulation
######################################################################

# load R/qtl and the map10 object
library(qtl)
data(map10)

# simulate under null
n.sim <- 10000
res0 <- rep(NA, n.sim)
for(i in 1:n.sim) {
  x <- sim.cross(map10[1:19], n.ind=250, type="f2")
  x <- calc.genoprob(x, step=1)
  out <- scanone(x, method="hk")
  res0[i] <- max(out[,3])
}

# estimated LOD threshold
print(thr <- quantile(res0, 0.95))

# genome length
print(G <- sum(summary(map10)[1:19,"length"]))

# Estimate significance threshold
thresh(G, "f2", d=10, p=0.05)

# histogram of Genome-wide maximum LOD scores
hist(res0, breaks=100, xlab="Genome-wide maximum LOD score")
rug(res0)

# simulate data with one QTL
alpha <- sqrt(2*0.08/(1-0.08))
n.sim <- 10000
loda <- est <- lo <- hi <- rep(NA, n.sim)
for(i in 1:n.sim) {
  x <- sim.cross(map10[1], n.ind=250, type="f2", 
                 model=c(1, 54, alpha, 0))
  x <- calc.genoprob(x, step=1)
  out <- scanone(x, method="hk")
  loda[i] <- max(out[,3])
  temp <- out[out[,3]==loda[i],2]
  if(length(temp) > 1) temp <- sample(temp, 1)
  est[i] <- temp
  li <- lodint(out)
  lo[i] <- li[1,2]
  hi[i] <- li[nrow(li),2]
}

# estimated power
mean(loda >= thr)

# Estimate power with R/qtlDesign
powercalc("f2", 250, sigma2=1, effect=c(alpha,0), thresh=thr,
          theta=0.09)

# histogram of genome-wide maximum LOD scores
hist(loda, breaks=100, xlab="Maximum LOD score")

# histogram of genome-wide maximum LOD scores
hist(est, breaks=100, xlab="Estimated QTL location (cM)")
rug(map10[[1]])

# Average estimated location
mean(est)

# SD of estimated locations
sd(est)

# SD of estimated location when LOD exceeds threshold
sig <- (loda >= thr)
sd(est[sig])

# Estimated coverage of 1.5-LOD support interval
mean(lo <= 54 & hi >= 54)

# Estimated coverage when LOD exceeds threshold
mean(lo[sig] <= 54 & hi[sig] >= 54)

# histogram of width of 1.5-LOD support interval
hist(hi[sig]-lo[sig], breaks=100, 
     xlab="Width of 1.5-LOD support interval (cM)")

# end of chap06.R
