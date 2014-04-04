######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap04.R"
#
# Chapter 4: Single-QTL analysis
#######################################################################

######################################################################
# Figure 4.1
#
#  Dot plots of the blood pressure phenotype against the
#  genotype at two selected markers, for the hyper data.
#  Confidence intervals for the average phenotype in each genotype
#  group are shown.
######################################################################
library(qtl)
data(hyper)
cx <- 0.6
set.seed(95109525)
par(mfrow=c(1,2), cex=cx, cex.lab=1/cx, cex.axis=1/cx)
par(mar=c(4.1,4.6,3.1,0.6))
plot.pxg(hyper, "D4Mit214")
par(mar=c(4.1,5.1,3.1,0.1))
plot.pxg(hyper, "D12Mit20")



######################################################################
# Figure 4.2
#
#  LOD scores for each marker on chromosomes 4 and 12 for the
#  hyper data, calculated by marker
#  regression.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.mr, chr=c(4, 12), ylab="LOD score")



######################################################################
# Figure 4.3
#
#  Phenotype distributions conditional on the genotype at two
#  markers, in the case of a backcross containing a single QTL.  The
#  markers are separated by 20 cM and the QTL sits within the marker
#  interval, 7 cM from the left marker. The dashed curves are the
#  distributions of the groups with common QTL genotype; the solid
#  curves are the mixture distributions.
######################################################################
par(mar=c(4.1,0.1,0.1,0.1))

x <- seq(0,100,by=0.1)
r <- (1-exp(-20/50))/2
r1 <- (1-exp(-7/50))/2
r2 <- (1-exp(-13/50))/2
p11 <- (1-r1)*(1-r2)/(1-r)
p12 <- (1-r1)*r2/r
p21 <- r1*(1-r2)/r
p22 <- r1*r2/(1-r)

par(bty="n")
plot(0,0,type="l",xlab="Phenotype",ylab="",xlim=c(20,100),ylim=c(0,5.3),yaxt="n")
a <- par("usr")
abline(h=c(0,1,2,3)*1.5 + a[3])
lines(x,(p11*dnorm(x,65,8)+(1-p11)*dnorm(x,50,8))/0.044+a[3])

text(rep(20,4),c(0,1,2,3)*1.5+a[3]+0.65, c("AB/AB","AB/AA","AA/AB","AA/AA"),
     adj=c(0,0.5))

lines(x,(p12*dnorm(x,65,8)+(1-p12)*dnorm(x,50,8))/0.044+1.5+a[3])
lines(x,p12*dnorm(x,65,8)/0.044+1.5+a[3],lty=2)
lines(x,(1-p12)*dnorm(x,50,8)/0.044+1.5+a[3],lty=2)

lines(x,(p21*dnorm(x,65,8)+(1-p21)*dnorm(x,50,8))/0.044+3+a[3])
lines(x,p21*dnorm(x,65,8)/0.044+3+a[3],lty=2)
lines(x,(1-p21)*dnorm(x,50,8)/0.044+3+a[3],lty=2)

lines(x,(p22*dnorm(x,65,8)+(1-p22)*dnorm(x,50,8))/0.044+4.5+a[3])

x <- seq(20,100,by=20)
segments(x,0.0+a[3],x,0.0+a[3]-0.15,xpd=TRUE, lend=1, ljoin=1)
segments(x,1.5+a[3],x,1.5+a[3]-0.15,xpd=TRUE, lend=1, ljoin=1)
segments(x,3.0+a[3],x,3.0+a[3]-0.15,xpd=TRUE, lend=1, ljoin=1)
segments(x,4.5+a[3],x,4.5+a[3]-0.15,xpd=TRUE, lend=1, ljoin=1)

z <- c(0,1.5,3,4.5)
segments(50,z+a[3],50,z+a[3]-0.15,lwd=2,xpd=TRUE, lend=1, ljoin=1, col="red")
segments(65,z+a[3],65,z+a[3]-0.15,lwd=2,xpd=TRUE, lend=1, ljoin=1, col="blue")
axis(side=1, at=51, labels=expression(mu[AA]), tick=FALSE, col.axis="red")
axis(side=1, at=66, labels=expression(mu[AB]), tick=FALSE, col.axis="blue")




######################################################################
# Figure 4.4
#
#  LOD scores by standard interval mapping for the hyper
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.em, ylab="LOD score")



######################################################################
# Figure 4.5
#
#  LOD scores for selected chromosomes for the hyper data
#  by standard interval mapping (blue) and marker regression
#  (red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.em, out.mr, chr=c(4, 12), col=c("blue", "red"),
     ylab="LOD score")



######################################################################
# Figure 4.6
#
#  LOD scores for selected chromosomes for the hyper data
#  by standard interval mapping (blue) and Haley--Knott regression
#  (red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.em, out.hk, chr=c(1,4,15), col=c("blue","red"),
     ylab="LOD score")



######################################################################
# Figure 4.7
#
#  Differences in the LOD scores from Haley--Knott regression and
#  standard interval mapping for selected chromosomes from the
#  hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.hk - out.em, chr=c(1,4,15), ylim=c(-0.5, 1.0),
     ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)



######################################################################
# Figure 4.8
#
#  LOD scores for selected chromosomes for the hyper data
#  by standard interval mapping (black), Haley--Knott regression (blue),
#  and the extended Haley--Knott method (red,
#  dashed).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.em, out.hk, out.ehk, chr=c(1,4,15), ylab="LOD score",
     lty=c(1,1,2))



######################################################################
# Figure 4.9
#
#  Differences in the LOD scores from Haley--Knott regression (in
#  blue) and the extended Haley--Knott method (in red) from the LOD
#  scores of standard interval mapping, for selected chromosomes from
#  the hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.hk - out.em, out.ehk - out.em, chr=c(1, 4, 15),
     col=c("blue", "red"), ylim=c(-0.5, 1),
     ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)



######################################################################
# Figure 4.10
#
#  Illustration of multiple imputations for a single backcross
#  individual.  Red and blue squares correspond to homozygous and
#  heterozygous genotypes, respectively, while open squares indicate
#  missing data.  The marker genotype data for the individual is shown
#  at the top, below a genetic map (in cM) for the chromosome; multiple
#  imputations of the genotype data, at the markers and at intervening
#  positions at 2 cM steps along the chromosome, are shown
#  below.
######################################################################

set.seed(12201969+15) 

data(hyper)
hyper <- subset(hyper, chr=12, ind=1:10)
m <- c(0,16,22,40,56)
names(m) <- names(hyper$geno[[1]]$map)
hyper$geno[[1]]$map <- m
marpos <- round(pull.map(hyper)[[1]])

par(mar=rep(0.1,4),bty="n")
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=c(6,100), ylim=c(3.5,103.5), xaxs="i", yaxs="i")
xl <- c(30,97)
yp <- 93
yd <- 1
segments(xl[1], yp, xl[2], yp)

L <- diff(range(marpos))
otherpos <- seq(0, L, by=2)
otherpos <- otherpos[is.na(match(otherpos, marpos))]

xotherpos <- otherpos * diff(xl)/L + xl[1]
xmarpos <- marpos * diff(xl)/L + xl[1]
segments(xmarpos, yp-yd*2, xmarpos, yp+yd*2)
segments(xotherpos, yp-yd, xotherpos, yp+yd)
text(xmarpos, yp+yd*6, marpos)

xd <- 4
text(xl[1]-xd, yp, "Genetic map:", adj=c(1,0.5))

yp2 <- yp-yd*8
ind <- 9 # or 13, 14, 57, 58
nmar <- length(marpos)
nother <- length(otherpos)
points(xmarpos, rep(yp2,nmar), pch=22, cex=1.5,
       bg=c("red", "blue")[hyper$geno[[1]]$data[ind,]], col="black")
points(xotherpos, rep(yp2, nother), pch=0, cex=1.5)

text(xl[1]-xd, yp2, "Observed data:", adj=c(1,0.5))

hyper <- sim.geno(hyper, n.draws=64, step=2)
hyper$geno[[1]]$draws <- hyper$geno[[1]]$draws[,,-(11:22)] 
xallpos <- sort(c(xmarpos, xotherpos))
npos <- length(xallpos)
yp3 <- yp2-yd*8
text(xl[1]-xd, yp3, "Imputations:", adj=c(1,0.5))
i <- 1
for(i in 1:15) {
  points(xallpos, rep(yp3,npos), pch=22, cex=1.5,
         bg=c("red", "blue")[hyper$geno[[1]]$draws[ind,,i]], col="black")
  yp3 <- yp3 - yd*5
}

yp4 <- yp2 - yd*6 - yd*3*5
xp2 <- 12
yd2 <- yd*6
points(xp2, yp4, pch=22, cex=1.5, bg="red", col="black")
text(xp2+xd*0.5, yp4, "= AA", adj=c(0,0.5))

points(xp2, yp4-yd2, pch=22, cex=1.5, bg="blue", col="black")
text(xp2+xd*0.5, yp4-yd2, "= AB", adj=c(0,0.5))

points(xp2, yp4-2*yd2, pch=0, cex=1.5)
text(xp2+xd*0.5, yp4-2*yd2, "= missing", adj=c(0,0.5))




######################################################################
# Figure 4.11
#
#  Illustration of the imputation results for chromosome 4 of
#  the \mboxhyper data.  The individual LOD curves from
#  n.draw imputations are shown in gray; the combined LOD score
#  from a total of n.imp imputations is shown in
#  black.
######################################################################

  set.seed(12201969)
  data(hyper)
  n.imp <- 64
  hyper <- subset(hyper,chr="4")
  hyper <- sim.geno(hyper, step=1, n.draws=n.imp)
  out.imp <- scanone(hyper, method="imp")
  n.draw <- 16
  temp <- as.data.frame(matrix(ncol=3+n.draw, nrow=nrow(out.imp)))
  temp[,1:3] <- out.imp
  names(temp)[1:3] <- names(out.imp)
  rownames(temp) <- rownames(out.imp)
  out.imp <- temp
  class(out.imp) <- c("scanone", "data.frame")
  temp <- hyper
  for(i in 1:n.draw) {
    temp$geno[[1]]$data <- hyper$geno[[1]]$draws[,,i]
    temp$geno[[1]]$map <- attr(hyper$geno[[1]]$draws,"map")
    out.imp[,i+3] <- scanone(temp, method="mr")[,3]
  }

par(mar=c(4.1,4.1,0.1,0.1))
plot(out.imp, lod=2, col="gray", ylab="LOD score",
     ylim=c(0,max(unlist(out.imp[,-(1:2)]))))
for(i in 2:n.draw)
  plot(out.imp, lod=i, add=TRUE, col="gray")
plot(out.imp, add=TRUE, col="black")




######################################################################
# Figure 4.12
#
#  LOD scores for selected chromosomes for the hyper data
#  by standard interval mapping (blue) and multiple imputation
#  (red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.em, out.imp, chr=c(1,4,15), col=c("blue", "red"),
     ylab="LOD score")



######################################################################
# Figure 4.13
#
#  Differences in the LOD scores from multiple imputation and
#  standard interval mapping for selected chromosomes from the
#  hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.imp - out.em, chr=c(1,4,15), ylim=c(-0.5, 0.5),
     ylab=expression(LOD[IMP] - LOD[EM]))
abline(h=0, lty=3)



######################################################################
# Figure 4.14
#
#  LOD scores by standard interval mapping (black), Haley--Knott
#  regression (blue) and the extended Haley--Knott method (red, dashed)
#  for chromosome 1 of the listeria data, when the genotype data
#  for all but the terminal markers have been
#  omitted.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
listeria <- calc.genoprob(listeria, step=1, error.prob=0.001)
outl.em <- scanone(listeria, chr=1)
outl.hk <- scanone(listeria, chr=1, method="hk")
outl.ehk <- scanone(listeria, chr=1, method="ehk")
plot(outl.em, outl.hk, outl.ehk, ylab="LOD score", 
     lty=c(1,1,2))



######################################################################
# Figure 4.15
#
#  Relationship between blood pressure and the genotype at
#  D4Mit214 for the hyper data.  The horizontal green segments
#  indicate the within-group averages.  In the left panel, the data for
#  all individuals is shown with the regression line of phenotype on
#  genotype.  In the right panel, the genotype data for all but the 92
#  individuals with extreme phenotypes is omitted.  The blue line is
#  the regression line obtained with only the 92 phenotyped
#  individuals.  The red line comes from a regression with all
#  individuals, but with those not genotyped placed intermediate
#  between the two genotype groups, as is done in Haley--Knott
#  regression.
######################################################################

cx <- 0.6
par(mfrow=c(1,2), cex=cx, cex.lab=1/cx, cex.axis=1/cx, las=1)
x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]
ex <- order(y)[c(1:46, 251-(1:45))]

par(mar=c(4.1,4.6,0.1,1.1))
plot(y ~ x, type="n", xlab="Genotype", ylab="bp",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n")
abline(v=1:2, lty=2, col="gray50")
points(x+u, y)
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="green2")
abline(lm(y~x)$coef, col="blue", lwd=2)

par(mar=c(4.1,5.6,0.1,0.1))
z <- x; z[-ex] <- 1.5
plot(y ~ z, type="n", xlab="Genotype", ylab="bp",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n")
abline(v=1:2, lty=2, col="gray50")
points(z+u, y, col=c("black", "gray70")[(z==1.5) + 1])
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, z, mean)
zz <- as.numeric(names(me))
segments(zz-0.15, me, zz+0.15, me, lwd=3, col="green2")
abline(lm(y[ex]~z[ex])$coef, col="blue", lwd=2)
abline(lm(y~z)$coef, col="red", lwd=2)




######################################################################
# Figure 4.16
#
#  LOD curves from standard interval mapping (black),
#  Haley--Knott regression (blue) and the extended Haley--Knott method
#  (red, dashed), for the hyper data when all individuals are
#  considered but genotype data for individuals with intermediate
#  phenotype are omitted. 
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out1.em, out1.hk, out1.ehk, ylab="LOD score", 
     lty=c(1,1,2))



######################################################################
# Figure 4.17
#
#  Differences in the LOD curves of Haley--Knott regression
#  (blue) and the extended Haley--Knott method (red), from the LOD
#  curves of standard interval mapping, for the hyper data when
#  all individuals are considered but genotype data for individuals
#  with intermediate phenotype are omitted. 
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out1.hk-out1.em, out1.ehk-out1.em, col=c("blue", "red"),
     ylim=c(-0.1, 4), ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)



######################################################################
# Figure 4.18
#
#  Differences in the LOD curves (for the hyper data)
#  from standard interval mapping (in black), Haley--Knott regression
#  (in green, dotted) and the extended Haley--Knott method (in red,
#  dashed), all calculated with only the data on the individuals with
#  extreme phenotypes, from the LOD curves of standard interval
#  mapping, with all individuals but with genotype data for individuals
#  with intermediate phenotype omitted.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out2.em-out1.em, out2.hk-out1.em, out2.ehk-out1.em,
     ylim=c(-0.1, 0.2), col=c("black", "green2", "red"), 
     lty=c(1,3,2), ylab="Difference in LOD scores")
abline(h=0, lty=3)



######################################################################
# Figure 4.19
#
#  Approximate null distribution of the LOD score at a
#  particular position (dashed curve) and of the maximum LOD score
#  genome-wide (solid curve) for a backcross in the
#  mouse.
######################################################################
  set.seed(44389856)
  x <- seq(0,5,length=513)
  y1 <- dchisq(x*2*log(10), 1)*(2*log(10))
  y2 <- density(apply(matrix(rchisq(15000000,1),ncol=150),1,max)/(2*log(10)),
                from=0, to=max(x), n=1024, width=0.4)
par(mar=c(4.1,0.1,0.1,0.1))
plot(x, y1, type="l", ylim=c(0,1), lwd=2, lty=2,
     xlab="LOD score", ylab="", yaxt="n")
lines(y2$x,y2$y,lwd=2)



######################################################################
# Figure 4.20
#
#  Diagram of the interval mapping process.
######################################################################

par(mar=rep(0.1,4), bty="n")
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(2,103),ylim=c(0,9), yaxs="i")
u <- par("usr")
x0 <- 5; x1 <- 25; y1 <- 8; xm <- mean(c(x0,x1))
segments(c(x0,x0,x0,x1),c(0,0,y1,0),c(x1,x0,x1,x1),c(0,y1,y1,y1))
text(xm,y1/2,"genotype\ndata")
y2 <- mean(c(u[4],y1))
text(xm,y2,"markers")
text(x0/2,y1/2,"individuals",srt=90)
x2 <- 27.5; x3 <- 32.5
segments(c(x2,x2,x2,x3),c(0,0,y1,0),c(x3,x2,x3,x3),c(0,y1,y1,y1))
text(mean(c(x2,x3)),y1/2,"phenotypes",srt=90)
x4 <- 50; x5 <- 70; y3 <- 3; y4 <- 5; xm <- mean(c(x4,x5)); ym <- mean(c(y3,y4))
segments(c(x4,x4,x4,x5),c(y3,y3,y4,y3),c(x5,x4,x5,x5),c(y3,y4,y4,y4))
text(xm,ym,"LOD scores")
arrows(x3+3, ym, x4-3, ym, length=0.1)
x6 <- 87.5; x7 <- 97.5
arrows(x5+3, ym, x6-3, ym, length=0.1)
text(x7,ym,"maximum\nLOD score")




######################################################################
# Figure 4.21
#
#  Histogram of the genome-wide maximum LOD scores from 1000
#  permutation replicates with the hyper data.  The LOD scores
#  were calculated by Haley--Knott regression.
######################################################################
par(mar=c(4.1,0.1,1.1,0.1), las=1)
hist(operm, ylab="", yaxt="n", breaks=seq(0, 5.0, len=101),
     xlab="maximum LOD score", main="")
rug(operm)



######################################################################
# Figure 4.22
#
#  The behavior of the X chromosome in a backcross.  Circles and
#  squares correspond to females and males, respectively.  Blue and
#  red bars correspond to DNA from strains A and B, respectively.  The
#  small bar is the Y chromosome.
######################################################################

lightblue <- "blue"
pink <- "red"

set.seed(122069)
par(mar=rep(0.1,4))
plot(0,0,type="n",xlim=c(2,72),ylim=c(105,0),xaxt="n",yaxt="n",
     xlab="",ylab="",xaxs="i",yaxs="i")
x <- 0; y <- 0

abline(v=37,h=52.5)

######################################################################
# (AxB)xA
######################################################################
cl <- 7
clx <- 2.5
cw <- 1.5*0.714
cd <- 0.5
cdown <- 4
xc <- x+9; xd <- 15
yc <- y+10; yd <- 12

text(19.5,y+4,"(A x B) x A",cex=1.4)
text(5,y+4,"a",cex=1.4,font=2)
segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1, col=lightblue, lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1, col=lightblue, lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1, col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1, col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1, col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1, col=lightblue,
     lend=1, ljoin=1)

ycp <- yc+yd
yd <- yd + 5

segments(xc+xd/2,yc,xc+xd/2,ycp,lwd=1)

xc <- xc
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"BC",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"BC",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1, col=lightblue,
     lend=1, ljoin=1)


######################################################################
# (BxA)xA
######################################################################
x <- 35; y <- 0
xc <- x+9
yc <- y+10;yd <- 12

text(x+19.5,y+4,"(B x A) x A",cex=1.4)
text(x+5,y+4,"b",cex=1.4,font=2)
segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1, col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

ycp <- yc+yd
yd <- yd + 5

segments(xc+xd/2,yc,xc+xd/2,ycp,lwd=1)

xc <- xc
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"BC",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"BC",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1, col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1, col=lightblue,
     lend=1, ljoin=1)


######################################################################
# Ax(AxB)
######################################################################
x <- 0; y <- 53
xc <- x+9 + xd/2
yc <- y+10;yd <- 12

text(19.5,y+4,"A x (A x B)",cex=1.4)
text(5,y+4,"c",cex=1.4,font=2)
segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc-xd/2
yc <- yc+yd

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[1]),cex=1.4)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

ycp <- yc+yd
yd <- yd + 5

segments(xc+xd/2,yc,xc+xd/2,ycp,lwd=1)

xc <- xc
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"BC",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"BC",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)


######################################################################
# Ax(BxA)
######################################################################
x <- 35; y <- 53
xc <- x+9 + xd/2
yc <- y+10;yd <- 12

text(x+19.5,y+4,"A x (B x A)",cex=1.4)
text(x+5,y+4,"d",cex=1.4,font=2)
segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc-xd/2
yc <- yc+yd

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[1]),cex=1.4)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

ycp <- yc+yd
yd <- yd + 5

segments(xc+xd/2,yc,xc+xd/2,ycp,lwd=1)

xc <- xc
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"BC",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"BC",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)




######################################################################
# Figure 4.23
#
#  The behavior of the X chromosome in an intercross.  Circles
#  and squares correspond to females and males, respectively.  Blue and
#  red bars correspond to DNA from strains A and B, respectively.  The
#  small bar is the Y chromosome.
######################################################################

lightblue <- "blue"
pink <- "red"

set.seed(41372)
par(mar=rep(0.1,4))
plot(0,0,type="n",xlim=c(0,98),ylim=c(105,0),xaxt="n",yaxt="n",
     xlab="",ylab="",xaxs="i",yaxs="i")
x <- 0; y <- 0

abline(v=49,h=52.5)

######################################################################
# (AxB)x(AxB)
######################################################################
cl <- 7
clx <- 2.5
cw <- 1.5
cd <- 0.5
cdown <- 4
xc <- x+7; xd <- 12
yc <- y+10; yd <- 14

text(25,y+4,"(A x B) x (A x B)",cex=1.4)
text(2.5,y+4,"a",cex=1.4,font=2)

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd*2

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc-xd*2,yc,lwd=1)
points(xc,yc,pch=15,cex=5,lwd=1,col="white")
points(xc,yc,pch=0,cex=5,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)

points(xc-xd*2,yc,pch=16,cex=6,lwd=1,col="white")
points(xc-xd*2,yc,pch=1,cex=6,lwd=1)
text(xc-xd*2,yc,expression(F[1]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc-xd*2-cd,yc+cdown,xc-xd*2-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc-xd*2+cd,yc+cdown,xc-xd*2+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

ycp <- yc+yd-2.5
yd <- yd + 1.5

segments(xc-xd,yc,xc-xd,ycp,lwd=1)

xc <- xc-xd*1.5
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[2]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[2]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)


######################################################################
# (BxA)x(AxB)
######################################################################
x <- 50; y <- 0
xc <- x+7; xd <- 12
yc <- y+10; yd <- 14

text(75,y+4,"(B x A) x (A x B)",cex=1.4)
text(52.5,y+4,"b",cex=1.4,font=2)

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd*2

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc-xd*2,yc,lwd=1)
points(xc,yc,pch=15,cex=5,lwd=1,col="white")
points(xc,yc,pch=0,cex=5,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)

points(xc-xd*2,yc,pch=16,cex=6,lwd=1,col="white")
points(xc-xd*2,yc,pch=1,cex=6,lwd=1)
text(xc-xd*2,yc,expression(F[1]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc-xd*2-cd,yc+cdown,xc-xd*2-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc-xd*2+cd,yc+cdown,xc-xd*2+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

ycp <- yc+yd-2.5
yd <- yd + 1.5

segments(xc-xd,yc,xc-xd,ycp,lwd=1)

xc <- xc-xd*1.5
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[2]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[2]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)


######################################################################
# (BxA)x(AxB)
######################################################################
x <- 0; y <- 53
xc <- x+7; xd <- 12
yc <- y+10; yd <- 14

text(25,y+4,"(A x B) x (B x A)",cex=1.4)
text(2.5,y+4,"c",cex=1.4,font=2)

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"A",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"B",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=pink,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd*2

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc-xd*2,yc,lwd=1)
points(xc,yc,pch=15,cex=5,lwd=1,col="white")
points(xc,yc,pch=0,cex=5,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)

points(xc-xd*2,yc,pch=16,cex=6,lwd=1,col="white")
points(xc-xd*2,yc,pch=1,cex=6,lwd=1)
text(xc-xd*2,yc,expression(F[1]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc-xd*2-cd,yc+cdown,xc-xd*2-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc-xd*2+cd,yc+cdown,xc-xd*2+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

ycp <- yc+yd-2.5
yd <- yd + 1.5

segments(xc-xd,yc,xc-xd,ycp,lwd=1)

xc <- xc-xd*1.5
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[2]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[2]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)


######################################################################
# (BxA)x(BxA)
######################################################################
x <- 50; y <- 53
xc <- x+7; xd <- 12
yc <- y+10; yd <- 14

text(75,y+4,"(B x A) x (B x A)",cex=1.4)
text(52.5,y+4,"d",cex=1.4,font=2)

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd*2

segments(xc,yc,xc+xd,yc,lwd=1)
points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,"B",cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,"A",cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd,yc+cdown,xc+xd-cd-cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

segments(xc+xd/2,yc,xc+xd/2,yc+yd,lwd=1)

xc <- xc+xd/2
yc <- yc+yd

segments(xc,yc,xc-xd*2,yc,lwd=1)
points(xc,yc,pch=15,cex=5,lwd=1,col="white")
points(xc,yc,pch=0,cex=5,lwd=1)
text(xc,yc,expression(F[1]),cex=1.4)

points(xc-xd*2,yc,pch=16,cex=6,lwd=1,col="white")
points(xc-xd*2,yc,pch=1,cex=6,lwd=1)
text(xc-xd*2,yc,expression(F[1]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)

rect(xc-xd*2-cd,yc+cdown,xc-xd*2-cd-cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc-xd*2+cd,yc+cdown,xc-xd*2+cd+cw,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)

ycp <- yc+yd-2.5
yd <- yd + 1.5

segments(xc-xd,yc,xc-xd,ycp,lwd=1)

xc <- xc-xd*1.5
yc <- yc+yd

segments(xc,ycp,xc+xd,ycp,lwd=1)
segments(xc,ycp,xc,yc,lwd=1)
segments(xc+xd,ycp,xc+xd,yc,lwd=1)

points(xc,yc,pch=16,cex=6,lwd=1,col="white")
points(xc,yc,pch=1,cex=6,lwd=1)
text(xc,yc,expression(F[2]),cex=1.4)
points(xc+xd,yc,pch=15,cex=5,lwd=1,col="white")
points(xc+xd,yc,pch=0,cex=5,lwd=1)
text(xc+xd,yc,expression(F[2]),cex=1.4)

rect(xc-cd-cw,yc+cdown,xc-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc-cd-cw,yc+cdown+cl*u1,xc-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+cd,yc+cdown,xc+cd+cw,yc+cdown+cl,lwd=1,col=pink,
     lend=1, ljoin=1)

rect(xc+xd-cd-cw,yc+cdown,xc+xd-cd,yc+cdown+cl,lwd=1,col=lightblue,
     lend=1, ljoin=1)
u <- runif(1,0.2,0.8)
u1 <- u; u2 <- 1
if(sample(1:2,1)==1) { u1 <- 0; u2 <- u }
rect(xc+xd-cd-cw,yc+cdown+cl*u1,xc+xd-cd,yc+cdown+cl*u2,lwd=1,col=pink,
     lend=1, ljoin=1)
rect(xc+xd+cd,yc+cdown,xc+xd+cd+cw,yc+cdown+clx,lwd=1,col=lightblue,
     lend=1, ljoin=1)




######################################################################
# Figure 4.25
#
#  Scatterplot of log_2(liver) versus log_2(spleen) for the
#  iron data, with females as red circles and males as blue
#  x's.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(log2(liver) ~ log2(spleen), data=iron$pheno, 
     col=c("red", "blue")[iron$pheno$sex],
     pch=c(1,4)[iron$pheno$sex])



######################################################################
# Figure 4.26
#
#  LOD curves for the liver phenotype in the iron data,
#  calculated by standard interval mapping.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(out.liver, ylab="LOD score")



######################################################################
# Figure 4.27
#
#  Illustration of the 1.5-LOD support interval for chromosome 4
#  of the hyper data.  The LOD curve was calculated by standard
#  interval mapping.
######################################################################


par(mar=c(4.1,4.1,0.1,0.1),xaxs="i",yaxs="i")
plot(out.em,ylim=c(0,8.5), chr=4, ylab="LOD score")
u <- par("usr")

li <- lodint(out.em,"4")
segments(li[2,2],li[2,3],40,li[2,3],lwd=2,lty=2,col="red")
segments(li[1,2],li[1,3],40,li[1,3],lwd=2,lty=2,col="red")
arrows(40,li[2,3]-0.2,40,li[1,3]+0.2,len=0.1,lwd=2,col="red")
text(42,mean(c(li[2,3],li[1,3])),"1.5", cex=0.8, adj=c(0,0.5), col="red")
u <- par("usr")
segments(li[1,2],u[3]+diff(u[3:4])*0.09,
         li[3,2],u[3]+diff(u[3:4])*0.09,lwd=2, col="blue")
for(i in c(1,3))
  segments(li[i,2],u[3]+diff(u[3:4])*0.07,
           li[i,2],u[3]+diff(u[3:4])*0.11,lwd=2, col="blue")
text(34,u[3]+diff(u[3:4])*0.09, "1.5-LOD support interval", adj=c(0,0.5),
     cex=0.8, col="blue")






######################################################################
# Figure 4.28
#
#  Illustration of the approximate 95% Bayes credible interval
#  for chromosome 4 of the hyper data.  The LOD curve was
#  calculated by standard interval mapping.
######################################################################


par(mar=c(4.1,4.6,0.1,0.1),xaxs="i", yaxs="i")

temp <- out.em[out.em[,1]=="4",]
loc <- temp[,2]
width <- diff((c(loc[1], loc) + c(loc, loc[length(loc)]))/2)
temp[,3] <- 10^temp[,3]
area <- temp[,3]*width
area <- area/sum(area)
temp[,3] <- area/width

plot(temp,ylim=c(0,0.9), chr=4, ylab=expression(paste(10^{LOD}, "   (rescaled)")))
u <- par("usr")

bi <- bayesint(out.em,"4")
left <- bi[1,2]
right <- bi[3,2]
wh <- temp[,2] >= left & temp[,2] <= right
x <- c(left, temp[wh,2], right)
y <- c(0, temp[wh,3], 0)
polygon(x, y, density=20, angle=45, col="red")
#segments(left, 0, left, temp[temp[,2]==left,3], lwd=2, col="red")
#segments(right, 0, right, temp[temp[,2]==right,3], lwd=2, col="red")
plot(temp, add=TRUE)
arrows(left, 0, left, 0.10, lwd=2, col="red", len=0.06)
arrows(right, 0, right, 0.10, lwd=2, col="red", len=0.06)
segments(left, 0, right, 0)

u <- par("usr")
segments(bi[1,2],u[3]+diff(u[3:4])*0.19,
         bi[3,2],u[3]+diff(u[3:4])*0.19,lwd=2, xpd=TRUE, col="blue")
for(i in c(1,3))
  segments(bi[i,2],u[3]+diff(u[3:4])*0.17,
           bi[i,2],u[3]+diff(u[3:4])*0.21,lwd=2, xpd=TRUE, col="blue")
text(34,u[3]+diff(u[3:4])*0.19, "95% Bayes interval", adj=c(0,0.5),
     cex=0.8, xpd=TRUE, col="blue")






######################################################################
# Figure 4.29
#
#  Histogram of the estimated QTL locations in 1000 bootstrap
#  replicates with the hyper data, chromosome 4, with LOD scores
#  calculated by standard interval mapping.  The tick marks beneath the
#  histogram indicate the locations of the genetic
#  markers.
######################################################################
par(mar=c(4.1,0.1,0.1,0.1),las=1)
hist(out.boot, breaks=brks, main="", xlab="Estimated QTL location",
     ylab="", yaxt="n")
rug(hyper$geno[[4]]$map)



######################################################################
# Figure 4.30
#
#  Illustration of selection bias in the estimated QTL effect.
#  The curves correspond to the distribution of the estimated percent
#  variance explained by a QTL for different values of the true effect
#  (indicated by the blue vertical lines).  The shaded regions
#  correspond to the cases where significant genomewide evidence for
#  the presence of a QTL would be obtained.  The red vertical lines
#  indicate the average estimated QTL effect, conditional on the
#  detection of the QTL.
######################################################################

  n.ind <- 250
  n.sim <- 100000

  hsq <- c(0.025, 0.05, 0.075)
  lod <- matrix(ncol=length(hsq), nrow=n.sim)
  a <- sqrt(hsq/(1-hsq))

  for(i in 1:n.sim) {
    for(j in seq(along=hsq)) {

      x <- sample(c(-1,1), n.ind, repl=TRUE)
      y <- x*a[j] + rnorm(n.ind)

      rss0 <- sum(lm(y ~ 1)$resid^2)
      rss1 <- sum(lm(y ~ factor(x))$resid^2)
      lod[i,j] <- (250/2)*log10(rss0/rss1)
    }
  }
  ve <- 1 - 10^(-2*lod/n.ind)

par(mfrow=c(3,1))
par(mar=c(5.1, 1.1, 3.1, 1.1))

for(i in 1:3) {
  temp <- hist(ve[,i], breaks=seq(0, 0.25, len=251), plot=FALSE)
  x <- rep(temp$breaks, rep(2, length(temp$breaks)))*100
  y <- c(0,rep(temp$density, rep(2, length(temp$counts))),0)
  truth <- c(2.5,5,7.5)[i]
  plot(x, y, main=paste("True variance explained = ", truth, "%", sep=""),
       yaxt="n", ylab="",       ylim=c(0,23), yaxs="i", type="l",
       xlab="Estimated percent variance explained", xlim=c(0,25), xaxs="i")

  thr <- 1-10^(-2/250*3)
  lo <- min(temp$breaks[temp$breaks >= thr])*100

  lod <- -250/2*log10(1-ve[,i])
  me <- mean(100*ve[lod>=3,i])
  power <- mean(lod>=3)*100
  bias <- (me - truth)/truth*100

  xx <- x[x>=lo]
  yy <- y[x>=lo]
  yy[1] <- 0
  polygon(xx,yy, density=10, lwd=2)
  abline(v=truth, lwd=2, col="blue")
  abline(v=me, lwd=2, col="red")

  text(20, 20, paste("Power = ", round(power), "%", sep=""), adj=c(0, 0.5), cex=1.3)
  text(20, 17, paste("Bias = ", round(bias), "%", sep=""), adj=c(0, 0.5), cex=1.3)
}



######################################################################
# Figure 4.31
#
#  Plot of the blood pressure against the genotype at marker
#  D4Mit164 for the hyper data.  The left panel was produced by
#  effectplot.  The right panel was produced by plot.pxg;
#  red dots correspond to imputed genotypes.  Error bars are +/- 1
#  SE.
######################################################################
set.seed(17748986)
cx <- 0.6
par(mfrow=c(1,2), cex=cx, cex.lab=1/cx, cex.axis=1/cx)
par(mar=c(5.1,4.6,0.1,1.1))
effectplot(hyper, mname1="D4Mit164", xlab="Genotype",
           main="")
par(mar=c(5.1,5.6,0.1,0.1))
plot.pxg(hyper, "D4Mit164", main="")



######################################################################
# Figure 4.32
#
#  LOD curves for selected chromosomes with the iron
#  data. Blue and red correspond to the liver and spleen phenotypes,
#  respectively.  Solid and dashed curves correspond to the original and
#  log scale, respectively.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
out.all <- scanone(iron, pheno.col=1:4)
plot(out.all, lodcolumn=1:2, col=c("blue", "red"), 
     chr=c(2, 7, 8, 9, 16), ylim=c(0,12.7), ylab="LOD score")
plot(out.all, lodcolumn=3:4, col=c("blue", "red"), lty=2,
     chr=c(2, 7, 8, 9, 16), add=TRUE)


# end of fig04.R
