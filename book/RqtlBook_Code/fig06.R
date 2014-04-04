######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap06.R"
#
# Chapter 6: Experimental design and power
#######################################################################

######################################################################
# Figure 6.1
#
#  Expected information from a selectively genotyped backcross
#  as a function of the selection fraction (the proportion of
#  extreme phenotypic individuals genotyped), in the middle of a marker
#  interval of length 0, 5, 10, and 20 cM.  The information is plotted
#  relative to a fully genotyped cross, where all individuals are
#  genotyped at a dense set of markers.  
######################################################################
x <- seq(0,1,by=0.01)
par(mar=c(3.6,4.1,2.1,2.6))
plot(x,info(x,0,"bc"),type="n",axes=F,xlim=c(0,1),ylim=c(0,1),
     xlab="",ylab="",yaxs="i",xaxs="i")
axis(1,at=seq(0,1,by=0.1),
     labels=as.character(seq(0,100,by=10)))
axis(2,at=seq(0,1,by=0.1),
     labels=as.character(seq(0,100,by=10)), las=1)
axis(3,at=seq(0,1,by=0.2),
     labels=as.character(seq(0,100,by=20)))
axis(4,at=seq(0,1,by=0.2),
     labels=as.character(seq(0,100,by=20)), las=1)
mtext(side=1,line=2,text="Selection fraction, in percent")
mtext(side=2,line=3,text="Expected percent of information")
lines(x,info(x,0,"bc"), lwd=2)
lines(c(0,0.4),c(info(0.4,0,"bc"),info(0.4,0,"bc")),lty=3)
lines(c(0,0.5),c(info(0.5,0,"bc"),info(0.5,0,"bc")),lty=3)
lines(c(0,0.6),c(info(0.6,0,"bc"),info(0.6,0,"bc")),lty=3)
lines(c(0.4,0.4),c(0,info(0.4,0,"bc")),lty=3)
lines(c(0.5,0.5),c(0,info(0.5,0,"bc")),lty=3)
lines(c(0.6,0.6),c(0,info(0.6,0,"bc")),lty=3)
lines(x,info(x,recomb(0.05),"bc"),lwd=2, col="blue")
lines(x,info(x,recomb(0.10),"bc"),lwd=2, col="red")
lines(x,info(x,recomb(0.20),"bc"),lwd=2, col="green2")
yloc <- c(info(0.9, recomb(0.05), "bc"),
          info(0.9, recomb(0.10), "bc"),
          info(0.9, recomb(0.20), "bc"))-0.02454
text(x=c(80,90,90,90)/100,y=c(0.972,yloc),
     labels=c("0 cM","5 cM","10 cM","20 cM"), 
     col=c("black", "blue", "red", "green2"))



######################################################################
# Figure 6.2
#
#  Illustration of selective phenotyping using simulated data.
#  We compare the LOD scores obtained using three different strategies:
#  the full cross of 200 individuals (black), 40 individuals
#  selectively phenotyped based on chromosome 1 genotypes (blue), and a
#  random set of 40 individuals (red).
######################################################################
mp <- sim.map(len=rep(100,5), n.mar=11, include.x=FALSE,
              eq.spacing=TRUE)
cr <- sim.cross(mp, model=c(1,50,1,0), n.ind=200, type="f2")
idx40 <- mma(pull.geno(cr,chr=1), p=40)
cr <- calc.genoprob(cr, step=2)
out1 <- scanone(cr)
out2 <- scanone(subset(cr,ind=idx40$cList))
out3 <- scanone(subset(cr,ind=sample(1:200,40)))
par(mar=c(4.1,4.1,0.1,0.1))
plot(out1, out2, out3, ylab="LOD score")



######################################################################
# Figure 6.3
#
#  Distribution of the genome-wide maximum LOD scores under the
#  null hypothesis of no QTL, for the case of an intercross with 250
#  individuals and with a genome modeled after that of the mouse and
#  with equally spaced markers at a 10 cM spacing.
######################################################################
par(mar=c(4.6,0.1,0.1,0.1))
hist(res0, breaks=seq(0, 8, len=100), 
     xlab="Genome-wide maximum LOD score",
     main="", yaxt="n", ylab="")
rug(res0)



######################################################################
# Figure 6.4
#
#  Distribution of the chromosome-wide maximum LOD scores in the
#  presence of a single QTL responsible for 8% of the phenotypic
#  variance, for the case of an intercross with 250 individuals and
#  with equally spaced markers at a 10 cM spacing.
######################################################################
par(mar=c(4.6,0.1,0.1,0.1))
hist(loda, breaks=seq(0, 15, len=100), 
     xlab="Maximum LOD score",
     main="", yaxt="n", ylab="")



######################################################################
# Figure 6.5
#
#  Estimated QTL location in 10,000 simulation replicates of an
#  intercross with 250 individuals, with a QTL located at 54 cM
#  (indicated by the blue triangle) and
#  responsible for 8% of the phenotypic variance. The tick marks at the bottom indicate the
#  marker locations (~10 cM spacing).
######################################################################
par(mar=c(4.6,0.1,0.1,0.1))
themap <- create.map(map10[[1]], 1, 0) # locations at which LOD scores calculated
br <- c(-0.5, (themap[-1] + themap[-length(themap)])/2,max(themap)+0.5) # breaks in histogram
hist(est, breaks=br,
     xlab="Estimated QTL location (cM)",
     main="", yaxt="n", ylab="", prob=TRUE)
rug(map10[[1]])
points(54, -0.0083, pch=17, col="blue", xpd=TRUE, cex=0.9)
text(54, -0.015, "QTL", col="blue", xpd=TRUE, cex=0.9)



######################################################################
# Figure 6.6
#
#  Distribution of the width of the 1.5-LOD support interval,
#  conditional on having significant evidence for a QTL, for the case
#  of a single QTL responsible for 8% of the phenotypic variance, in
#  an intercross with 250 individuals and markers at a 10 cM
#  spacing.
######################################################################
par(mar=c(4.6,0.1,0.1,0.1))
hist(hi[sig]-lo[sig], breaks=seq(0, 127, len=128),
     xlab="Width of 1.5-LOD support interval (cM)",
     main="", yaxt="n", ylab="")


# end of fig06.R
