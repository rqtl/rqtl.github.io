######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 1: Introduction
#######################################################################

######################################################################
# Figure 1.1
#
#  Schematic representation of the autosomes in a backcross
#  experiment.  The two inbred strains, A and B, are represented by
#  blue and pink chromosomes, respectively.  The F_1 generation,
#  obtained by crossing the two strains, receives a single chromosome
#  from each parent, and all individuals are genetically identical.  If
#  we cross an individual from the F_1 generation ``back'' to one of
#  the parental strains (the A strain, in this example), we obtain a
#  population exhibiting genetic variation.  The backcross individuals
#  receive an intact A chromosome from their A parent.  The chromosome
#  received from their F_1 parent may be an intact A or B chromosome,
#  but is generally a mosaic of the A and B chromosomes as a result of
#  recombination at meiosis.  Any given locus has a 50% chance of
#  being heterozygous and a 50% chance of being homozygous.  This
#  figure represents the autosomes only.  When considering the X
#  chromosome, four backcross populations (``directions'') are possible
#  (see Fig. ??? on
#  page ???).
######################################################################
# load utility functions
source("meiosis_func.R")

lightblue <- "blue"
pink <- "red"

par(mar=rep(0.1,4),bty="n")
plot(0,0,xlim=c(-15,864),ylim=c(0,480),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
xd <- 2
rect(c(300-xd+25,328-xd+25),c(480,480),c(310+xd+25,338+xd+25),c(385,385), lend=1, ljoin=1,
     col=lightblue)
rect(c(526-xd-25,554-xd-25),c(480,480),c(536+xd-25,564+xd-25),c(385,385),
     col=pink, lend=1, ljoin=1)

points(432,440,pch=4,cex=1.5)
#segments(432,400,432,340)
#segments(432,400,432,340)
#segments(319,340,545,340)
#arrows(c(319,545),c(340,340),c(319,545),c(300,300),len=0.1)
arrows(432, 400, 432, 300, len=0.1)

text(300,(480+385)/2,"A",cex=1.5,adj=c(1,0.5))
text(564,(480+385)/2,"B",cex=1.5,adj=c(0,0.5))

rect(413-xd,287,423+xd,192, lend=1, ljoin=1, col=lightblue)
rect(441-xd,287,451+xd,192, col=pink, lend=1, ljoin=1)

rect(413-176-xd,287,423-176+xd,192, lend=1, ljoin=1, col=lightblue)
rect(441-176-xd,287,451-176+xd,192, lend=1, ljoin=1, col=lightblue)


points(344,247,pch=4,cex=1.5)
segments(344,207,344,147)
segments(57,147,849,147)
arrows(seq(57,849,by=88),rep(147,10),seq(57,849,by=88),rep(107,10),len=0.1)

text(451+25,(287+192)/2,expression(F[1]),cex=1.5,adj=c(0,0.5))
text(413-176-25,(287+192)/2,"A",cex=1.5,adj=c(1,0.5))

f1 <- create.par(100,c(1,2))
set.seed(99019)
f2 <- vector("list",10)
for(i in 1:10) f2[[i]] <- cross(f1,f1,m=10,obl=TRUE)

xloc <- 38
mult <- 95/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
for(i in 1:10) {
  rect(xloc-xd,0,xloc+10+xd,95, lend=1, ljoin=1, col=lightblue)
  rect(xloc+28-xd,0,xloc+38+xd,95, lend=1, ljoin=1, col=lightblue)

#  f2m <- f2[[i]]$mat
#  for(j in 2:ncol(f2m)) {
#    if(f2m[2,j]==2)
#      rect(xloc-xd,f2m[1,j]*mult,xloc+10+xd,f2m[1,j-1]*mult,
#           density=20,angle=45, lend=1, ljoin=1)
#  }
  f2p <- f2[[i]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+28-xd,f2p[1,j]*mult,xloc+38+xd,f2p[1,j-1]*mult,
           col=pink, lend=1, ljoin=1)
  }
  xloc <- xloc+38+50
}
text(38-25,95/2,"BC",cex=1.5,adj=c(1,0.5))




######################################################################
# Figure 1.2
#
#  Schematic representation of the autosomes in an intercross
#  experiment.  The two inbred strains, A and B, are represented by
#  blue and pink chromosomes, respectively.  As with the backcross
#  strategy, the F_1 generation, obtained by crossing the two
#  strains, has a chromosome from each parent, and all F_1
#  individuals are genetically identical.  By crossing two individuals
#  in the F_1 generation (or by selfing when possible), we create
#  genetic variation in the resulting F_2 population.  The three
#  possible genotypes AA, AB, and BB appear in a 1:2:1 ratio.  This
#  figure represents the autosomes only; the behavior of the X
#  chromosome is displayed in Fig. ??? on
#  page ???.
######################################################################
# load utility functions
source("meiosis_func.R")

lightblue <- "blue"
pink <- "red"

par(mar=rep(0.1,4), bty="n")
plot(0,0,xlim=c(-15,864),ylim=c(0,480),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
xd <- 2
rect(c(300-xd,328-xd)+25,c(480,480),c(310+xd,338+xd)+25,c(385,385), lend=1, ljoin=1,
     col=lightblue)
rect(c(526-xd,554-xd)-25,c(480,480),c(536+xd,564+xd)-25,c(385,385),
     col=pink, lend=1, ljoin=1)

points(432,440,pch=4,cex=1.5)
segments(432,400,432,340)
segments(319+25,340,545-25,340)
arrows(c(319+25,545-25),c(340,340),c(319+25,545-25),c(300,300),len=0.1)

text(300,(480+385)/2,"A",cex=1.5,adj=c(1,0.5))
text(564,(480+385)/2,"B",cex=1.5,adj=c(0,0.5))

rect(300-xd+25,287,310+xd+25,192, lend=1, ljoin=1, col=lightblue)
rect(328-xd+25,287,338+xd+25,192, col=pink, lend=1, ljoin=1)
rect(526-xd-25,287,536+xd-25,192, lend=1, ljoin=1, col=lightblue)
rect(554-xd-25,287,564+xd-25,192, col=pink, lend=1, ljoin=1)

points(432,247,pch=4,cex=1.5)
segments(432,207,432,147)
segments(57,147,849,147)
arrows(seq(57,849,by=88),rep(147,10),seq(57,849,by=88),rep(107,10),len=0.1)

text(300,(287+192)/2,expression(F[1]),cex=1.5,adj=c(1,0.5))
text(564,(287+192)/2,expression(F[1]),cex=1.5,adj=c(0,0.5))

f1 <- create.par(100,c(1,2))
set.seed(999)
f2 <- vector("list",10)
for(i in 1:10) f2[[i]] <- cross(f1,f1,m=10,obl=TRUE)

xloc <- 38
mult <- 95/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
for(i in 1:10) {
  rect(xloc-xd,0,xloc+10+xd,95, lend=1, ljoin=1, col=lightblue)
  rect(xloc+28-xd,0,xloc+38+xd,95, lend=1, ljoin=1, col=lightblue)

  f2m <- f2[[i]]$mat
  for(j in 2:ncol(f2m)) {
    if(f2m[2,j]==2)
      rect(xloc-xd,f2m[1,j]*mult,xloc+10+xd,f2m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f2p <- f2[[i]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+28-xd,f2p[1,j]*mult,xloc+38+xd,f2p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  xloc <- xloc+38+50
}
text(38-25,95/2,expression(F[2]),cex=1.5,adj=c(1,0.5))




######################################################################
# Figure 1.3
#
#  Schematic representation of the autosomes during the breeding
#  of recombinant inbred lines by sibling mating.  The first two
#  generations are identical to an F_2 intercross.  In subsequent
#  generations, siblings are mated producing progeny that are less and
#  less heterozygous.  If continued indefinitely, this process will
#  produce individuals that are completely homozygous at every locus,
#  but with chromosomes that are a mosaic of the parental chromosomes.
#  The frequency of the breakpoints between the AA and BB genotypes is
#  determined by the breeding scheme (sib-mating or selfing, or some
#  other scheme). In practice, 10--20 generations of inbreeding are
#  actually performed.
######################################################################
# load utility functions
source("meiosis_func.R")

lightblue <- "blue"
pink <- "red"

par(mar=rep(0.1,4), bty="n")
plot(0,0,xlim=c(-40,864),ylim=c(30,480),xaxt="n",yaxt="n",xlab="",ylab="",type="n")
u <- par("usr")

# initial generation
xd <- 1
rect(c(300-xd,328-xd)+25,c(480,480),c(310+xd,338+xd)+25,c(440,440), lend=1, ljoin=1,
     col=lightblue)
rect(c(526-xd,554-xd)-25,c(480,480),c(536+xd,564+xd)-25,c(440,440),
      col=pink,lend=1, ljoin=1)
points(432,465,pch=4,cex=1.5)
segments(432,450,432,430)
segments((310+328)/2+25,430,(554+536)/2-25,430)
arrows(c((310+328)/2+25,(554+536)/2-25),c(430,430),
       c((310+328)/2+25,(554+536)/2-25),c(410,410),len=0.1)
text(300-xd+25-32,460,"A",adj=c(1,0.5),cex=1.5)
text(564+xd-25+32,460,"B",adj=c(0,0.5),cex=1.5)

# F1 hybrids
rect(300-xd+25,400,310+xd+25,360, lend=1, ljoin=1, col=lightblue)
rect(526-xd-25,400,536+xd-25,360, lend=1, ljoin=1, col=lightblue)
rect(328-xd+25,400,338+xd+25,360, col=pink, lend=1, ljoin=1)
rect(554-xd-25,400,564+xd-25,360, col=pink, lend=1, ljoin=1)
text(300-xd+25-32,380,expression(F[1]),adj=c(1,0.5),cex=1.5)
points(432,385,pch=4,cex=1.5)
segments(432,370,432,350)


# simulate the recombinant inbred lines (by sib mating)
  f1 <- create.par(100,c(1,2))
  set.seed(994)
  f2 <- vector("list",10)
  for(i in 1:10) f2[[i]] <- cross(f1,f1,m=10,obl=TRUE)
  f3 <- vector("list",10)
  for(i in 1:5) {
    f3[[2*i-1]] <- cross(f2[[2*i-1]],f2[[2*i]],m=10,obl=TRUE)
    f3[[2*i]] <- cross(f2[[2*i-1]],f2[[2*i]],m=10,obl=TRUE)
  }
  f4 <- vector("list",10)
  for(i in 1:5) {
    f4[[2*i-1]] <- cross(f3[[2*i-1]],f3[[2*i]],m=10,obl=TRUE)
    f4[[2*i]] <- cross(f3[[2*i-1]],f3[[2*i]],m=10,obl=TRUE)
  }
  temp <- f4
  my.ri8 <- vector("list",10)
  for(j in 1:60) {
    for(i in 1:5) {
      my.ri8[[2*i-1]] <- cross(temp[[2*i-1]],temp[[2*i]],m=10,obl=TRUE)
      my.ri8[[2*i]] <- cross(temp[[2*i-1]],temp[[2*i]],m=10,obl=TRUE)
    }
    temp <- my.ri8
  }

# F2 generation
xloc <- 38
mult <- 40/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
xxloc <- NULL
for(i in 1:4) {
  rect(xloc-6,320,xloc+6,280, lend=1, ljoin=1, col=lightblue)
  rect(xloc+20,320,xloc+32,280, lend=1, ljoin=1, col=lightblue)
  rect(xloc+76,320,xloc+88,280, lend=1, ljoin=1, col=lightblue)
  rect(xloc+102,320,xloc+114,280, lend=1, ljoin=1, col=lightblue)

  f2m <- f2[[2*i-1]]$mat
  for(j in 2:ncol(f2m)) {
    if(f2m[2,j]==2)
      rect(xloc-6,280+f2m[1,j]*mult,xloc+6,280+f2m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f2p <- f2[[2*i-1]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+20,280+f2p[1,j]*mult,xloc+32,280+f2p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  f2m <- f2[[2*i]]$mat
  for(j in 2:ncol(f2m)) {
    if(f2m[2,j]==2)
      rect(xloc+76,280+f2m[1,j]*mult,xloc+88,280+f2m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f2p <- f2[[2*i]]$pat
  for(j in 2:ncol(f2p)) {
    if(f2p[2,j]==2)
      rect(xloc+102,280+f2p[1,j]*mult,xloc+114,280+f2p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  xxloc <- c(xxloc,xloc+(6+20)/2,xloc+(88+102)/2)

  points(xloc+54,300,pch=4,cex=1.5)
  segments(xloc+54,290,xloc+54,270)
  segments(xxloc[2*i-1],270,xxloc[2*i],270)
  
  xloc <- xloc+78+120+30
}
segments(min(xxloc),350,max(xxloc),350)
arrows(xxloc,c(350,350),xxloc,c(330,330),len=0.1)
text(0,300,expression(F[2]),adj=c(1,0.5),cex=1.5)
arrows(xxloc,270,xxloc,250,len=0.1)

# F3 generation
xloc <- 38
for(i in 1:4) {
  rect(xloc-6,240,xloc+6,200, lend=1, ljoin=1, col=lightblue)
  rect(xloc+20,240,xloc+32,200, lend=1, ljoin=1, col=lightblue)
  rect(xloc+76,240,xloc+88,200, lend=1, ljoin=1, col=lightblue)
  rect(xloc+102,240,xloc+114,200, lend=1, ljoin=1, col=lightblue)

  f3m <- f3[[2*i-1]]$mat
  for(j in 2:ncol(f3m)) {
    if(f3m[2,j]==2)
      rect(xloc-6,200+f3m[1,j]*mult,xloc+6,200+f3m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f3p <- f3[[2*i-1]]$pat
  for(j in 2:ncol(f3p)) {
    if(f3p[2,j]==2)
      rect(xloc+20,200+f3p[1,j]*mult,xloc+32,200+f3p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  f3m <- f3[[2*i]]$mat
  for(j in 2:ncol(f3m)) {
    if(f3m[2,j]==2)
      rect(xloc+76,200+f3m[1,j]*mult,xloc+88,200+f3m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f3p <- f3[[2*i]]$pat
  for(j in 2:ncol(f3p)) {
    if(f3p[2,j]==2)
      rect(xloc+102,200+f3p[1,j]*mult,xloc+114,200+f3p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  points(xloc+54,220,pch=4,cex=1.5)
  segments(xloc+54,210,xloc+54,190)
  segments(xxloc[2*i-1],190,xxloc[2*i],190)
  
  xloc <- xloc+78+120+30
}
text(0,220,expression(F[3]),adj=c(1,0.5),cex=1.5)
arrows(xxloc,190,xxloc,170,len=0.1)

# F4 generation
xloc <- 38
for(i in 1:4) {
  rect(xloc-6,160,xloc+6,120, lend=1, ljoin=1, col=lightblue)
  rect(xloc+20,160,xloc+32,120, lend=1, ljoin=1, col=lightblue)
  rect(xloc+76,160,xloc+88,120, lend=1, ljoin=1, col=lightblue)
  rect(xloc+102,160,xloc+114,120, lend=1, ljoin=1, col=lightblue)

  f4m <- f4[[2*i-1]]$mat
  for(j in 2:ncol(f4m)) {
    if(f4m[2,j]==2)
      rect(xloc-6,120+f4m[1,j]*mult,xloc+6,120+f4m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f4p <- f4[[2*i-1]]$pat
  for(j in 2:ncol(f4p)) {
    if(f4p[2,j]==2)
      rect(xloc+20,120+f4p[1,j]*mult,xloc+32,120+f4p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  f4m <- f4[[2*i]]$mat
  for(j in 2:ncol(f4m)) {
    if(f4m[2,j]==2)
      rect(xloc+76,120+f4m[1,j]*mult,xloc+88,120+f4m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  f4p <- f4[[2*i]]$pat
  for(j in 2:ncol(f4p)) {
    if(f4p[2,j]==2)
      rect(xloc+102,120+f4p[1,j]*mult,xloc+114,120+f4p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  points(xloc+54,140,pch=4,cex=1.5)
  arrows(xloc+54,87,xloc+54,80,len=0.1)
  arrows(xloc+54,130,xloc+54,80,len=0.1,lty=3)
  
  xloc <- xloc+78+120+30
}
text(0,140,expression(F[4]),adj=c(1,0.5),cex=1.5)

# final RIL
xloc <- 38
a <- 70-40
mult <- (70-a)/f2[[1]]$mat[1,ncol(f2[[1]]$mat)]
for(i in 1:4) {
  rect(xloc+54-19,70,xloc+54-7,a, lend=1, ljoin=1, col=lightblue)
  rect(xloc+54+7,70,xloc+54+19,a, lend=1, ljoin=1, col=lightblue)

  my.ri8m <- my.ri8[[2*i-1]]$mat
  for(j in 2:ncol(my.ri8m)) {
    if(my.ri8m[2,j]==2)
      rect(xloc+54-19,a+my.ri8m[1,j]*mult,xloc+54-7,a+my.ri8m[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }
  my.ri8p <- my.ri8[[2*i-1]]$pat
  for(j in 2:ncol(my.ri8p)) {
    if(my.ri8p[2,j]==2)
      rect(xloc+54+7,a+my.ri8p[1,j]*mult,xloc+54+19,a+my.ri8p[1,j-1]*mult,
            col=pink, lend=1, ljoin=1)
  }

  xloc <- xloc+78+120+30
}
points(rep(-16,3),c(-8,0,8)+mean(c(130,(70+a)/2)),pch=16,cex=0.2)
text(0,(70+a)/2,expression(F[infinity]),adj=c(1,0.5),cex=1.5)



######################################################################
# Figure 1.4
#
#  (Hypothetical) phenotype data for the two parental strains,
#  the F_1 hybrid, and a backcross population.  Vertical line
#  segments are plotted at the within-group averages.  In this
#  simulated example, the difference in the phenotype means between
#  strains is quite large relative to the within-strain variation.  The
#  mean phenotype in the F_1 generation is intermediate between the
#  two strains, while the backcross progeny exhibit a wider spectrum of
#  phenotypic variation and resemble the A strain more than the B
#  strain.
######################################################################
library(qtl)
set.seed(4131972)
data(map10)
bc <- sim.cross(map10, type="bc", n.ind=50,
                model=cbind(1:5, 50, 0.5))
bc$pheno <- ((2*bc$qtlgeno - 3) %*% rep(1, 5)) + rnorm(50, 31.7*0.75+45*0.25, 2.5)

x <- rev(list("A" = rnorm(15, 31.7, 2.5),
          "B" = rnorm(15, 45.4, 2.5),
          "F1" = rnorm(15, (31.7+45.4)/2,  2.5),
          "BC" = bc$pheno[,1]))
par(las=1, mar=c(4.1,2.1,0.1,0.1))
stripchart(x, method="jitter", pch=1, xlab="Phenotype")
abline(h=1:4, lty=2)
me <- sapply(x, mean)
segments(me, 1:4-0.2, me, 1:4+0.2, lwd=3, col="blue")



######################################################################
# Figure 1.5
#
#  Histogram of systolic blood pressure for 250 backcross mice
#  from Sugiyama et al. (2001).  Also shown are the phenotypic ranges of parental
#  and F_1 hybrid strains (mean +/- 2 SD).  The A/J (or A) strain
#  is normotensive, while the C57BL/6J (or B6) strain is hypertensive.
#  The F_1 hybrids and the backcross resemble the hypertensive B6
#  strain. Notice that the range of the phenotype is not unlike
#  humans.
######################################################################
library(qtl) 
data(hyper)
par(mar=c(4.1,0.1,0.1,0.1))
hist(hyper$pheno[,1], breaks=seq(70,130,by=2), las=1, main="",
     xlab="Blood pressure (mm of Hg)",ylim=c(0,40),
     axes=FALSE, ylab="")
axis(1,at=seq(80,130,by=10))
strain.smry <- c(104.8,7.4)
hashmark <- function(m,s,pos,adj=0,label,  ...)
{
  lines(c(m-2*s,m+2*s),c(pos,pos), ...)
  points(m,pos, pch=16, ...)
  text(m+2*s+1,pos+adj,label,pos=4, ...)
}
hashmark(104.8,7.4,39,0,"C57BL/6J", col="red")
hashmark(85.0,5.9,36,0,"A/J", col="blue")
hashmark(101.1,7.5,33,0.1,expression(paste(group("(",B6xA,")"),F[1])), col="green2")
hashmark(100.6,5.0,30,-1.0,expression(paste(group("(",AxB6,")"),F[1])), col="orange3")



######################################################################
# Figure 1.6
#
#  The genetic map of markers typed in the data from
#  Sugiyama et al. (2001).  Almost all marker intervals are less than 20cM.  Some
#  regions, most notably on chromosomes 1 and 4, have a higher density
#  of markers.
######################################################################
hyper <- drop.nullmarkers(hyper)
hyper <- jittermap(hyper)
par(mar=c(4.1,4.1,0.1,0.1))
plot.map(hyper, main="", alternate.chrid=TRUE)



######################################################################
# Figure 1.8
#
#  The statistical structure of the QTL mapping problem.  The
#  QTL and covariates are responsible for phenotypic variation
#  (indicated by the directed solid arrows).  The markers and the QTL
#  are correlated with each other due to linkage (indicated by the
#  bidirectional solid arrow).  The markers do not directly cause the
#  phenotype; some markers may be associated with the phenotype via
#  linkage to the QTL (indicated by the directed dashed
#  arrow).
######################################################################
par(mar=rep(0.1,4), cex=1.3, bty="n")
plot(0,0, type="n", xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=c(12,88), ylim=c(9,41))
text(20,10,"Markers")
text(60,10,"Phenotype")
text(40,40,"QTL")
text(80,40,"Covariates")
arrows(28, 10, 50, 10, lty=2, len=0.1, lwd=2)
arrows(49, 10, 50, 10, len=0.1, lwd=2)
arrows(22, 15, 38, 35, len=0.1, code=3, lwd=2)
arrows(42, 35, 58, 15, len=0.1, lwd=2)
arrows(78, 35, 62, 15, len=0.1, lwd=2)



######################################################################
# Figure 1.9
#
#  Illustration of the problem of inferring missing genotype
#  data.  Each row is the marker genotype data for a different
#  backcross individual.  Dashes indicate missing data.  We seek the
#  probability that each individual has genotype AA or AB at the
#  putative QTL indicated by the triangle.
######################################################################
par(mar=rep(0.1,4), bty="n")
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",
     xlim=c(-3,103),ylim=c(-2.9,11.1), yaxs="i", xaxs="i")
x <- c(0,19.7,32.8,64.2,85.4,100)
segments(0,9,100,9)
segments(x,8.5,x,9.5)
points(55,9.6,pch=25, col="orange3", bg="orange3")
text(x, 10.5, c(expression(M[1]), expression(M[2]), expression(M[3]),
                expression(M[4]), expression(M[5]), expression(M[6])))
text(55, 10.6, expression(Q), col="orange3")
text(x,6.5,c("AA","AA","AA","AA","AB","AB"), 
     col=c("blue","blue","blue","blue","red","red"))
text(x,4,c("AB","AA","AA","-","AB","AA"),
     col=c("red","blue","blue","black","red","blue"))
text(x,1.5,c("AA","AB","-","AB","AA","AA"),
     col=c("blue","red","black","red","blue","blue"))
arrows(55,-1.5,55,0, len=0.1, col="orange3")
text(55,-2.5,"?", col="orange3")



######################################################################
# Figure 1.10
#
#  Illustration of possible effects of two QTL in a backcross.
#  A. Additive QTL.  B. Interacting QTL.  Points are
#  located at the average phenotype for a given two-locus genotype.
#  Line segments connect the averages for a given genotype at the
#  second QTL.
######################################################################

me1 <- c(10,40,  60,90)
g1 <- c(1,2,1,2)
g2 <- g1 + 3
me2 <- c(10,40,  25,90)

par(mar=c(3.6,4.1,0.1,0.1))
plot(g1, me1, lwd=2, xaxt="n", xlab="", las=1,
     xlim=c(0.5,6.5), ylim=c(0,100), ylab="Ave. phenotype",
     xaxs="i")
abline(v=3.5)
axis(side=1, at=1:2, labels=c("AA", "AB"))
lines(g1[1:2], me1[1:2], lwd=2, col="blue")
lines(g1[3:4], me1[3:4], lwd=2, col="red")
points(g1, me1, col="white", pch=16, lwd=2)
points(g1, me1, lwd=2, col=c("blue","blue","red","red"))
u <- par("usr")
text(1.5, u[3]-diff(u[3:4])*0.23, "QTL 1", xpd=TRUE, cex=1.1)
text(2.3, me1[c(2,4)], c("AA", "AB"), col=c("blue","red"))
text(3, me1[4], "QTL 2", cex=1.1)
text(0.7, 97, "A", font=2, cex=1.2)
text(3.7, 97, "B", font=2, cex=1.2)

points(g2, me2, lwd=2)
axis(side=1, at=1:2+3, labels=c("AA", "AB"))
lines(g2[1:2], me2[1:2], lwd=2, col="blue")
lines(g2[3:4], me2[3:4], lwd=2, col="red")
points(g2, me2, col="white", pch=16, lwd=2)
points(g2, me2, lwd=2, col=c("blue","blue","red","red"))
text(4.5, u[3]-diff(u[3:4])*0.23, "QTL 1", xpd=TRUE, cex=1.1)
text(5.3, me2[c(2,4)], c("AA", "AB"), col=c("blue","red"))
text(6, me2[4], "QTL 2", cex=1.1, xpd=TRUE)




######################################################################
# Figure 1.11
#
#  Illustration of possible effects of two QTL in an intercross.
#  A. Additive QTL.  B. Interacting QTL.  Points are
#  located at the average phenotype for a given two-locus genotype.
#  Line segments connect the averages for a given genotype at the
#  second QTL.
######################################################################

me1 <- c(10,20,60,  20,30,70,  40,50,90)
g1 <- rep(c(1,1.5,2), 3)
g2 <- g1 + 3
me2 <- c(10,20, 0,  20,30,30,  40,50,90)

par(mar=c(3.6,4.1,0.1,0.1))
plot(g1, me1, lwd=2, xaxt="n", xlab="", las=1,
     xlim=c(0.5,6.5), ylim=c(0,100), ylab="Ave. phenotype",
     xaxs="i")
abline(v=3.5)
axis(side=1, at=g1[1:3], labels=c("AA", "AB", "BB"))
lines(g1[1:3], me1[1:3], lwd=2, col="blue")
lines(g1[4:6], me1[4:6], lwd=2, col="red")
lines(g1[7:9], me1[7:9], lwd=2, col="green2")
points(g1, me1, col="white", pch=16, lwd=2)
points(g1, me1, lwd=2, col=rep(c("blue","red","green2"), rep(3,3)))
u <- par("usr")
text(1.5, u[3]-diff(u[3:4])*0.23, "QTL 1", xpd=TRUE, cex=1.1)
text(2.3, me1[c(3,6,9)], c("AA", "AB", "BB"), col=c("blue","red","green2"))
text(3, me1[9], "QTL 2", cex=1.1)
text(0.7, 97, "A", font=2, cex=1.2)
text(3.7, 97, "B", font=2, cex=1.2)

points(g2, me2, lwd=2)
axis(side=1, at=g2[1:3], labels=c("AA", "AB", "BB"))
lines(g2[1:3], me2[1:3], lwd=2, col="blue")
lines(g2[4:6], me2[4:6], lwd=2, col="red")
lines(g2[7:9], me2[7:9], lwd=2, col="green2")
points(g2, me2, col="white", pch=16, lwd=2)
points(g2, me2, lwd=2, col=rep(c("blue","red","green2"), rep(3,3)))
u <- par("usr")
text(4.5, u[3]-diff(u[3:4])*0.23, "QTL 1", xpd=TRUE, cex=1.1)
text(5.3, me2[c(3,6,9)], c("AA", "AB", "BB"), col=c("blue","red","green2"))
text(6, me2[9], "QTL 2", cex=1.1)



# end of fig01.R
