######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap10.R"
#
# Chapter 10: Case study I
#######################################################################

######################################################################
# Figure 10.2
#
#  Scatterplot of on1 and on2, the gonad-specific
#  ovariole counts in the ovar data, with points randomly
#  jittered so that overlapping points may be distinguished.
######################################################################
set.seed(68996116)
par(mar=c(4.1,4.1,0.1,0.1), las=1, pty="s")
plot(jitter(on2) ~ jitter(on1), data=ovar$pheno,
     xlab="on1", ylab="on2", cex=0.6,
     xlim=c(6.66, 17.53), ylim=c(6.66, 17.53))



######################################################################
# Figure 10.3
#
#  Box plot of onm for the two crosses forming the
#  ovar data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
boxplot(onm ~ cross, data=ovar$pheno, horizontal=TRUE,
        xlab="Ovariole count", ylab="Cross")



######################################################################
# Figure 10.4
#
#  Plot of ovariole count against number of typed marker
#  genotypes for the initial cross in the ovar
#  data.  On the right are points corresponding to the 94 genotyped
#  individuals; on the left are the remaining individuals, for which
#  only five morphological markers were scored.
######################################################################
set.seed(58184589)
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(jitter(ntyped(ovar1)),jitter(pull.pheno(ovar1,"onm")),
     xlab="No. genotypes", ylab="Ovariole count")



######################################################################
# Figure 10.5
#
#  Plot of estimated recombination fractions (upper left) and
#  LOD scores for a test of r=1/2 (lower right) for all pairs of
#  markers in the initial cross in the ovar data.  Red indicates
#  linkage, blue indicates no linkage, and gray indicates missing
#  values (for marker pairs that were not both typed in any one
#  individual).
######################################################################
par(mar=c(4.1,4.1,1.1,0.1), las=1, pty="s")
ovar1 <- est.rf(ovar1)
plot.rf(ovar1, main="")



######################################################################
# Figure 10.6
#
#  The genetic map in the ovar data plotted against the
#  map estimated from the individuals in the initial cross.  For each
#  chromosome, the line on the left is the map provided with the data,
#  and the line on the right is the map estimated using the Haldane mapping function; line
#  segments connect the positions for each marker.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
  newmap <- est.map(ovar1, error.prob=0.001)
plot.map(ovar1, newmap, main="")



######################################################################
# Figure 10.8
#
#  LOD curves from a genome scan by multiple imputation for the
#  initial cross in the ovar data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(out1, ylab="LOD score")



######################################################################
# Figure 10.9
#
#  Estimated phenotype averages for the two groups defined by
#  genotypes at marker cpo, for the initial cross in the
#  ovar1 data.  Error bars are +/- 1 SE.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, pty="s")
effectplot(ovar1, mname1="cpo", main="")



######################################################################
# Figure 10.10
#
#  LOD curves from a genome scan by multiple imputation for the
#  initial cross in the ovar data, adjusting for an inferred QTL
#  at the marker cpo.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(out1.c3, ylab="LOD score")



######################################################################
# Figure 10.11
#
#  Profile LOD curves for a three-QTL model, for the initial
#  cross in the ovar data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex=0.8)
plotLodProfile(rqtl, col=c("black","red","blue"), 
               ylab="Profile LOD score", xlim=c(-12.5, 270), 
               xaxs="i")



######################################################################
# Figure 10.12
#
#  Plot of the onm phenotype against genotype at the two
#  tightly linked loci on chr 2, for the initial cross in the
#  ovar data.  Red dots in the right panel are for imputed
#  genotypes.  Error bars are +/- 1 SE.
######################################################################
mar <- find.marker(ovar1, 2, c(92.4, 94))
par(mfrow=c(1,2))
par(mar=c(4.1,4.1,0.1,1.1))
effectplot(ovar1, mname1="2@92.4", mname2="2@94",
           ylim=c(9.5,16.5), main="", add.legend=FALSE)
u <- par("usr")
legend(c(u[1]+diff(u[1:2])*0.05,
         u[1]+diff(u[1:2])*0.45),
       c(16.2, 14.7), c("II", "IE"), lty=1, pch=1, cex=1,
       col=c("red","blue"), xjust=0.5, yjust=0.5)
       
text(u[1]+diff(u[1:2])*0.25, 16.5, "Amy-d")
par(mar=c(4.1,5.1,0.1,0.1))
set.seed(30022337)
plot.pxg(ovar1, marker=mar,  col=c("blue","green2"),
           ylim=c(9.5,16.5), main="", cex=0.6)
axis(side=1, 0.1, paste(mar[1], ":", sep=""), line=-0.5, xpd=TRUE, tick=FALSE)
axis(side=1, 0.1, paste(mar[2], ":", sep=""), line=0.5, xpd=TRUE, tick=FALSE)



######################################################################
# Figure 10.13
#
#  LOD curves from a genome scan by multiple imputation for the
#  ovar data.  The black curves are for the combined data, the
#  blue curves are for the initial cross, and the red curves are for
#  the second cross.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(outc, out1, out2, ylab="LOD score")



######################################################################
# Figure 10.14
#
#  LOD curves from a genome scan by multiple imputation for the
#  ovar data, adjusting for an inferred QTL on chr 3.  The black
#  curves are for the combined data, the blue curves are for the
#  initial cross, and the red curves are for the second cross. Note
#  that the position of the inferred QTL (on which we are conditioning)
#  is different for the initial cross versus for the combined data and
#  for the second cross.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
plot(outc.c3, out1.c3, out2.c3, ylab="LOD score")



######################################################################
# Figure 10.15
#
#  LOD scores from a two-dimensional, two-QTL scan on chr 3,
#  controlling for a locus on chr 2, for the ovar data.
#  lod_i is displayed in the upper left triangle; lod_fv1 is
#  displayed in the lower right triangle.  In the color scale on the
#  right, numbers to the left and right correspond to lod_i and
#  lod_fv1, respectively.
######################################################################
par(cex.axis=0.9)
plot(out.ap.c3, lower="cond-int", 
     mar1=c(4.1,4.1,0.6,1.6), mar2=c(4.1,2.1,0.6,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 10.16
#
#  LOD scores from a two-dimensional, two-QTL scan on chr 2,
#  controlling for a locus on chr 3, for the ovar data.
#  lod_av1 is displayed in the upper left triangle; lod_fv1
#  is displayed in the lower right triangle.  In the color scale on the
#  right, numbers to the left and right correspond to lod_av1 and
#  lod_fv1, respectively.
######################################################################
par(cex.axis=0.9)
plot(out.ap.c2, lower="cond-int",  upper="cond-add",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 10.17
#
#  Profile LOD curves for a four-QTL model, for the full
#  ovar data; the two QTL on chr 3 were allowed to
#  interact.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex=0.8)
plotLodProfile(rqtl, col=c("red","blue","red","blue"),
               ylab="Profile LOD score")


# end of fig10.R
