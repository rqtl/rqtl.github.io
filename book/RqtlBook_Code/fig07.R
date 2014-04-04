######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap07.R"
#
# Chapter 7: Working with covariates
#######################################################################

######################################################################
# Figure 7.1
#
#  Illustration of the effects of a QTL and an additive
#  covariate in a backcross in the case of (A) sex as the covariate and
#  (B) a quantitative covariate.
######################################################################

me1 <- c(10,30,50,70)
g1 <- c(1,2,1,2)
g2 <- g1 + 3
me2 <- c(10,40,  25,90)

par(mar=c(3.6,4.1,0.1,0.1))
plot(g1, me1, lwd=2, xaxt="n", xlab="", las=1,
     xlim=c(0.5,4.5), ylim=c(0,80), ylab="Average phenotype",
     xaxs="i")
abline(v=2.5)
axis(side=1, at=1:2, labels=c("female", "male"))
lines(g1[1:2], me1[1:2], lwd=2)
lines(g1[3:4], me1[3:4], lwd=2, lty=2)
points(g1, me1, col="white", pch=16, lwd=2)
points(g1, me1, lwd=2)
u <- par("usr")
text(1.5, u[3]-diff(u[3:4])*0.23, "sex", xpd=TRUE, cex=1.1)
text(2.3, me1[c(2,4)], c("AA", "AB"))
text(0.7, 77, "A", font=2, cex=1.2)
text(2.7, 77, "B", font=2, cex=1.2)

arrows(1, me1[1]+2, 1, me1[3]-2, col="blue", len=0.1, code=3)
text(1.1, mean(me1[c(1,3)]), expression(beta[g]), col="blue")

arrows(2, me1[1]+1, 2, me1[2]-2, col="red", len=0.1, code=3)
segments(1.95, me1[1], 2.05, me1[1], col="red")
text(2.1, mean(me1[c(1,2)]), expression(beta[x]), col="red")

u <- par("usr")

text(4.25, 40, "AA")
text(4.25, 80, "AB")
thex <- seq(2.7, 4.3, len=5)
axis(side=1, at=thex, label=0:4)
text((2.5+u[2])/2, u[3]-diff(u[3:4])*0.23, "x", xpd=TRUE, cex=1.1)
arrows(2.7, (2.7-3)*20+10 + 1, 2.7, (2.7-3)*20+50 - 1, col="blue",
       len=0.1, code=3)
text(2.8, mean(c((2.7-3)*20+10, (2.7-3)*20+50)), expression(beta[g]), col="blue")
segments(thex[3], (thex[3]-3)*20+10, thex[4], (thex[3]-3)*20+10, col="red")
segments(thex[4], (thex[3]-3)*20+10, thex[4], (thex[4]-3)*20+10, col="red")
text(thex[4]+0.1, ((thex[3]-3)*20+10 + (thex[4]-3)*20+10)/2, expression(beta[x]), col="red")

segments(2.5, (2.5-3)*20+10, u[2], (u[2]-3)*20+10, lwd=2)
segments(2.5, (2.5-3)*20+50, u[2], (u[2]-3)*20+50, lwd=2, lty=2)






######################################################################
# Figure 7.2
#
#  Illustration of the effect of a QTL as a function of w, in
#  the model implied by the use of y/w as the phenotype in QTL
#  mapping.
######################################################################

par(mar=c(4.1,10.1,0.1,6.1))
plot(0, 0, xlab="w", las=1, type="n", yaxs="i",
     xlim=c(0,4.5), ylim=c(0,0.35), ylab="Average phenotype",
     xaxs="i", yaxt="n")
axis(side=2, at=seq(0, 0.3, by=0.1), las=1)
segments(2, 2*0.03, 3, 2*0.03, col="blue")
segments(3, 2*0.03, 3, 3*0.03, col="blue")
text(3.1, 2.5*0.03, expression(mu), col="blue", adj=c(0,0.5))

segments(2, 2*0.07, 3, 2*0.07, col="red")
segments(3, 2*0.07, 3, 3*0.07, col="red")
text(3.1, 2.5*0.07, expression(mu+beta[g]), col="red", adj=c(0,0.5))

text(4, 0.14, "AA")
text(4, 0.31, "AB")

abline(0, 0.03, lwd=2)
abline(0, 0.07, lty=2, lwd=2)





######################################################################
# Figure 7.3
#
#  Diagram of the interval mapping process in the presence of
#  additive covariates.
######################################################################

par(mar=rep(0.1,4), bty="n")
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(2,113),ylim=c(0,9), yaxs="i")
u <- par("usr")
x0 <- 5; x1 <- 25; y1 <- 8; xm <- mean(c(x0,x1))
segments(c(x0,x0,x0,x1),c(0,0,y1,0),c(x1,x0,x1,x1),c(0,y1,y1,y1))
text(xm,y1/2,"genotype\ndata")
y2 <- mean(c(u[4],y1))
text(xm,y2,"markers")
text(x0/2,y1/2,"individuals",srt=90)
x2 <- 27.5; x3 <- 32.5; x3b <- x3+(x3-x2)*2
segments(c(x2,x2,x2,x3),c(0,0,y1,0),c(x3,x2,x3,x3),c(0,y1,y1,y1))
text(mean(c(x2,x3)),y1/2,"phenotypes",srt=90)
segments(c(x3b,x3b,x3b,x3),c(0,0,y1,0),c(x3,x3b,x3,x3),c(0,y1,y1,y1))
text(mean(c(x3,x3b)),y1/2,"covariates",srt=90)
x4 <- x3b+17.5; x5 <- x4+20; y3 <- 3; y4 <- 5; xm <- mean(c(x4,x5)); ym <- mean(c(y3,y4))
segments(c(x4,x4,x4,x5),c(y3,y3,y4,y3),c(x5,x4,x5,x5),c(y3,y4,y4,y4))
text(xm,ym,"LOD scores")
arrows(x3b+3, ym, x4-3, ym, length=0.1)
x6 <- x5+17.5; x7 <- x6+10
arrows(x5+3, ym, x6-3, ym, length=0.1)
text(x7,ym,"maximum\nLOD score")




######################################################################
# Figure 7.4
#
#  Box plots of the gut length phenotype by sex and cross in the
#  gutlength data.
######################################################################
par(las=1, mar=c(4.1, 10.6, 0.1, 1.1))
boxplot(gutlength ~ sex*cross, data=gutlength$pheno, 
        horizontal=TRUE, xlab="Gut length (cm)",
        col=c("red","blue"))



######################################################################
# Figure 7.5
#
#  Plot of LOD scores for gut length with no covariates (in
#  blue) and with inclusion of sex and cross as additive covariates (in
#  red, dashed) for the gutlength data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.0, out.a, col=c("blue", "red"), lty=1:2, 
     ylab="LOD score", alternate.chrid=TRUE)



######################################################################
# Figure 7.6
#
#  Plot of differences between the LOD scores for gut length
#  with sex and cross as additive covariates versus without the
#  covariates for the gutlength data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.a - out.0, ylab="LOD w/ covar - LOD w/o covar", 
     ylim=c(-1, 1), alternate.chrid=TRUE)
abline(h=0, lty=2)



######################################################################
# Figure 7.7
#
#  Illustration of the effects of a QTL and an interactive
#  covariate in a backcross in the case of (A) sex as the covariate and
#  (B) a quantitative covariate.
######################################################################

thegreen <- "green2"

me1 <- c(30,20,50,70)
g1 <- c(1,2,1,2)
g2 <- g1 + 3
me2 <- c(10,40,  25,90)

par(mar=c(3.6,4.1,0.1,0.1))
plot(g1, me1, lwd=2, xaxt="n", xlab="", las=1,
     xlim=c(0.5,4.5), ylim=c(10,80), ylab="Average phenotype",
     xaxs="i")
abline(v=2.5)
axis(side=1, at=1:2, labels=c("female", "male"))
lines(g1[1:2], me1[1:2], lwd=2)
lines(g1[3:4], me1[3:4], lwd=2, lty=2)
points(g1, me1, col="white", pch=16, lwd=2)
points(g1, me1, lwd=2)
u <- par("usr")
text(1.5, u[3]-diff(u[3:4])*0.23, "sex", xpd=TRUE, cex=1.1)
text(2.3, me1[c(2,4)], c("AA", "AB"))
text(0.7, 77, "A", font=2, cex=1.2)
text(2.7, 77, "B", font=2, cex=1.2)

arrows(1, me1[1]+2, 1, me1[3]-2, col="blue", len=0.1, code=3)
text(1.05, mean(me1[c(1,3)]), expression(beta[g]), col="blue",
     adj=c(0, 0.5))

arrows(1, me1[1]-2, 1, me1[2]+1, col="red", len=0.1, code=3)
segments(0.95, me1[2], 1.05, me1[2], col="red")
text(1.05, mean(me1[c(1,2)]), expression(beta[x]), col="red",
     adj=c(0, 0.5))

arrows(2, me1[2]+2, 2, me1[4]-2, col=thegreen, len=0.1, code=3)
text(2.05, mean(me1[c(2,4)]), expression(beta[g]+gamma), col=thegreen,
     adj=c(0, 0.5))


u <- par("usr")

text(4.25, 22.5, "AA")
text(4.25, 80, "AB")
thex <- seq(2.7, 4.3, len=5)
axis(side=1, at=thex, label=0:4)
text((2.5+u[2])/2, u[3]-diff(u[3:4])*0.23, "x", xpd=TRUE, cex=1.1)
arrows(2.7, -(2.7-3)*10+30 + 1, 2.7, (2.7-3)*20+50 - 1, col="blue",
       len=0.1, code=3)
text(2.75, mean(c(-(2.7-3)*10+30, (2.7-3)*20+50)), expression(beta[g]), col="blue",
     adj=c(0, 0.5))

segments(thex[3], -(thex[3]-3)*10+30, thex[4], -(thex[3]-3)*10+30, col="red")
segments(thex[4], -(thex[3]-3)*10+30, thex[4], -(thex[4]-3)*10+30, col="red")
text(thex[4]+0.05, (-(thex[3]-3)*10+30 - (thex[4]-3)*10+30)/2, expression(beta[x]), col="red",
     adj=c(0, 0.5))

segments(thex[3], (thex[3]-3)*20+50, thex[4], (thex[3]-3)*20+50, col=thegreen)
segments(thex[4], (thex[3]-3)*20+50, thex[4], (thex[4]-3)*20+50, col=thegreen)
text(thex[4]+0.05, ((thex[3]-3)*20+50 + (thex[4]-3)*20+50)/2, expression(beta[x]+gamma), col=thegreen,
     adj=c(0, 0.5))

segments(2.5, -(2.5-3)*10+30, u[2], -(u[2]-3)*10+30, lwd=2)
segments(2.5, (2.5-3)*20+50, u[2], (u[2]-3)*20+50, lwd=2, lty=2)





######################################################################
# Figure 7.8
#
#  Plot of lod_f (in blue), with cross as an additive
#  covariate and sex as an interactive covariate, and lod_i (in
#  red), for the QTL x sex interaction, for the gutlength
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.i, out.i - out.a, ylab="LOD score", 
     col=c("blue", "red"), alternate.chrid=TRUE)



######################################################################
# Figure 7.9
#
#  Plot of LOD scores for males (in blue) and females (in red)
#  for the gutlength data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.m, out.f, col=c("blue", "red"), ylab="LOD score", 
     alternate.chrid=TRUE)



######################################################################
# Figure 7.10
#
#  Plot of the difference between the sum of the LOD scores from
#  the separate analyses of males and females and lod_f, from their
#  joint analysis, for the gutlength data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.m + out.f - out.i, ylim=c(-0.5,0.5), 
     ylab="LOD(males) + LOD(females) - LODf", 
     alternate.chrid=TRUE)
abline(h=0, lty=2)



######################################################################
# Figure 7.11
#
#  Plot of the estimated phenotype averages +/- 1 SE as a
#  function of sex (with males in blue and females in red) and
#  genotype, at the positions nearest the peak LOD score on four selected
#  chromosomes, for the gutlength data. C and B correspond to
#  the C3HeBFeJ and C57BL/6J alleles, respectively.
######################################################################
gutlength <- sim.geno(subset(gutlength,chr=c(4,5,18,"X")),
                      n.draws=128, step=1, error.prob=0.001)
par(mfrow=c(2,2), mar=c(4.1,4.1,1.6,0.6))
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="4@65", main="Chromosome 4", 
           add.legend=FALSE)
par(mar=c(4.1,4.6,1.6,0.1))
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="5@22", main="Chromosome 5", 
           add.legend=FALSE)
par(mar=c(4.1,4.1,1.6,0.6))
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="18@48.2", main="Chromosome 18", 
           add.legend=FALSE)
par(mar=c(4.1,4.6,1.6,0.1))
effectplot(gutlength, mname1="sex", ylim=c(15.1, 17.2),
           mname2="X@58", main="X chromosome", 
           add.legend=FALSE)



######################################################################
# Figure 7.12
#
#  Plot of lod_f (in black), lod_a (in blue), and lod_i
#  (in red), for the nf1 data, with parent-of-origin of the
#  NPcis mutation considered as a covariate.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.all, lod=1:3, ylab="LOD score")



######################################################################
# Figure 7.13
#
#  LOD scores for the analysis of the nf1 data, split by
#  parent-of-origin of the NPcis mutation, with results for
#  individuals receiving the mutation from their mother and father in
#  red and blue, respectively.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.bypoo, lod=1:2, col=c("red", "blue"), 
     ylab="LOD score")



######################################################################
# Figure 7.14
#
#  Proportion of affecteds as a function of genotype and the
#  parent-of-origin of the NPcis mutation, for the nf1
#  data.
######################################################################
nf1 <- sim.geno(subset(nf1,chr=c(15,19)), n.draws=128, 
                step=1, error.prob=0.001)
par(mfrow=c(1,2), yaxs="i")
par(mar=c(4.1,4.1,0.6,0.6))
effectplot(nf1, mname1="NPcis", mark1=1-from.mom, main="",
           geno1=c("Mom", "Dad"), ylim=c(0,1), mname2="15@13",
           legend.lab="")
par(mar=c(4.1,4.6,0.6,0.1))
effectplot(nf1, mname1="NPcis", mark1=1-from.mom, main="",
           geno1=c("Mom", "Dad"), ylim=c(0,1), mname2="19@0",
           legend.lab="")



######################################################################
# Figure 7.15
#
#  LOD scores for the analysis of the hyper data, with no
#  covariates (in blue) and with imputed genotypes at D4Mit164 included
#  as an additive covariate (in red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out, out.ag, col=c("blue", "red"), ylab="LOD score")



######################################################################
# Figure 7.16
#
#  Interaction LOD scores, indicating evidence for an
#  interaction between a QTL and the marker D4Mit164, for the
#  hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.ig - out.ag, ylab="interaction LOD score")



######################################################################
# Figure 7.17
#
#  Results of composite interval mapping (CIM) for the
#  hyper data, with forward selection to three markers and three
#  different choices of window sizes.  In each panel, the results from
#  standard interval mapping are in blue, and those from CIM are in
#  red.  Selected marker covariates are indicated by green dots.
######################################################################
par(mar=c(4.1,4.1,2.1,0.1))
chr <- c(1, 2, 4, 6, 15)
par(mfrow=c(3,1))
plot(out, out.cim.20, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = 20 cM")
add.cim.covar(out.cim.20, chr=chr, bg="green2", pch=21, cex=1.5, col="black")
plot(out, out.cim.40, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = 40 cM")
add.cim.covar(out.cim.40, chr=chr, bg="green2", pch=21, cex=1.5, col="black")
plot(out, out.cim.inf, chr=chr, ylab="LOD score",
     col=c("blue", "red"), main="window = Inf")
add.cim.covar(out.cim.inf, chr=chr, bg="green2", pch=21, cex=1.5, col="black")


# end of fig07.R
