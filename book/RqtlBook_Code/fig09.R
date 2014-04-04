######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap09.R"
#
# Chapter 9: Fit and exploration of multiple-QTL models
#######################################################################

######################################################################
# Figure 9.1
#
#  An example of a regression tree.  This is a type of decision
#  tree that describes a QTL model.  Individuals with genotype AB or BB
#  at QTL 1 have average phenotype 20, individuals with genotype AA at
#  both QTL 1 and QTL 2 have average phenotype 100, and individuals
#  with genotype AA at QTL 1 and either AB or BB at QTL 2 have average
#  phenotype 80.
######################################################################
par(mar=c(0,0,0,0))
plot(0,0,type="n", xlab="", xaxt="n", ylab="", yaxt="n", bty="n",
     xlim=c(-10,95), ylim=c(5,96))
segments(c(10,30,50,50),c(10,50,90,90),c(30,50,30,70),c(50,10,50,50),lwd=2)
points(c(10,30,50,50,70),c(10,50,90,10,50),cex=6,lwd=2,pch=c(15,16,16,15,15),col="white")
points(c(10,30,50,50,70),c(10,50,90,10,50),cex=6,lwd=2,pch=c(15,16,16,15,15)-15)
text(c(10,50,70),c(10,10,50),c("100","80","20"), cex=1.3)
text(50,90,expression(q[1]),cex=1.3)
text(30,50,expression(q[2]),cex=1.3)
text(19,35,"AA",adj=1)
text(41,35,"AB or BB",adj=0)
text(39,75,"AA",adj=1)
text(61,75,"AB or BB",adj=0)



######################################################################
# Figure 9.2
#
#  Graphical depictions of QTL models, with nodes
#  corresponding to QTL and edges indicating interactions.  A:
#  Two interacting QTL.  B: Three QTL, of which two interact.
#  C: Three QTL with all pairs interacting.  D: A
#  seven-QTL model, with one QTL exhibiting pairwise interactions with
#  each of the other QTL.
######################################################################
par(mar=rep(0,4))
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n",
     xlim=c(10,110), ylim=c(3.5,92.5), xaxs="i", yaxs="i")

x <- c(25,40)
y <- c(75,85)
points(x, y, pch=16, cex=1.5)
segments(x[1], y[1], x[2], y[2], lwd=2, ljoin=1, lend=1)
text(45, 90, "A", font=2)

x <- c(20, 65/2, 38)
y <- c(15, 35, 10)
points(x, y, pch=16, cex=1.5)
segments(x[1], y[1], x[3], y[3], lwd=2, ljoin=1, lend=1)
text(45, 40, "B", font=2)

x <- c(20, 65/2, 38)+50
y <- c(15, 35, 10)+50
points(x, y, pch=16, cex=1.5)
segments(x[1], y[1], x[2], y[2], lwd=2, ljoin=1, lend=1)
segments(x[1], y[1], x[3], y[3], lwd=2, ljoin=1, lend=1)
segments(x[2], y[2], x[3], y[3], lwd=2, ljoin=1, lend=1)
text(95, 90, "C", font=2)

x <- c(80, 70, 66, 71, 87, 93, 91)-4
y <- c(20, 5, 23, 36, 5, 17, 33)
points(x, y, pch=16, cex=1.5)
for(i in 2:length(x))
  segments(x[1], y[1], x[i], y[i], lwd=2, ljoin=1, lend=1)
text(95, 40, "D", font=2)



######################################################################
# Figure 9.3
#
#  Locations of the QTL object qtl on the genetic map for
#  the hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), cex=0.7, cex.axis=1.2, cex.lab=1.2)
plot(qtl, main="")



######################################################################
# Figure 9.4
#
#  LOD profiles for a four-QTL model with the hyper
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), cex=0.8)
plotLodProfile(rqtl, ylab="Profile LOD score")



######################################################################
# Figure 9.5
#
#  LOD curves for adding one QTL to the four-QTL model, with the
#  hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.aq, ylab="LOD score")



######################################################################
# Figure 9.7
#
#  Interaction LOD curves in the scan for an additional QTL,
#  interacting with the chr 15 locus, to be added to the four-QTL
#  model, with the hyper data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.aqi - out.aq, ylab="LOD score")



######################################################################
# Figure 9.8
#
#  Results of a two-dimensional, two-QTL scan on chr 1, in the
#  context of a model with additional QTL on chr 4, 6 and 15, and a
#  6x15 interaction, with the hyper data. lod_av1 is
#  in the upper left triangle, and lod_fv1 is in the lower right
#  triangle.  In the color scale on the right, numbers to the left and
#  right correspond to lod_av1 and lod_fv1,
#  respectively.
######################################################################
plot(out.ap, lower="cond-int", upper="cond-add", 
     layout=list(cbind(1,2),c(4.5,1)),
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1))



######################################################################
# Figure 9.9
#
#  Results of a two-dimensional, two-QTL scan on chr 7 and 15,
#  in the context of a model with additional QTL on chr 1, 4, and 6,
#  with the hyper data.  The two QTL being scanned were allowed
#  to interact, and the first of them interacts with the chr 6 locus.
#  The LOD scores displayed are for the five-QTL model relative to the
#  three-QTL model.  The x-axis corresponds to the first of the new QTL
#  (which interacts with the chr 6 locus); the y-axis corresponds to
#  the second of the new QTL.
######################################################################
plot(out.ap2, 
     layout=list(cbind(1,2),c(4.5,1)),
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1))



######################################################################
# Figure 9.10
#
#  The sequence of models visited by the forward/backward search
#  of stepwiseqtl, with the hyper
#  data.
######################################################################
par(mar=c(0.6,0.1,2.1,0.6))
par(mfrow=c(6,3))
for(i in seq(along=thetrace))
  plotModel(thetrace[[i]], chronly=TRUE,
            main=paste(i, ": pLOD =", 
              round(attr(thetrace[[i]], "pLOD"), 2)))


# end of fig09.R
