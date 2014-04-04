######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap11.R"
#
# Chapter 11: Case study II
#######################################################################

######################################################################
# Figure 11.2
#
#  Box plots of the tth phenotype (time to hatch)
#  according to the female source of the egg for the individuals, for
#  the trout data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1)
boxplot(tth ~ female, data=trout$pheno, 
        xlab="Female", ylab="Time to hatch")



######################################################################
# Figure 11.3
#
#  Plot of estimated recombination fractions (upper left) and
#  LOD scores for a test of r=1/2 (lower right) for all pairs of
#  markers for the trout data.  Red indicates linkage; blue
#  indicates no linkage.
######################################################################
par(mar=c(4.1,4.1,2.1,2.1), las=1, pty="s", cex.axis=0.8)
trout <- est.rf(trout)
plot.rf(trout, alternate.chrid=TRUE, main="")



######################################################################
# Figure 11.4
#
#  Plot of estimated recombination fractions (upper left) and
#  LOD scores for a test of r=1/2 (lower right) for all pairs of
#  markers on selected linkage groups, for the trout data.  Red
#  indicates linkage; blue indicates no linkage.
######################################################################
par(mar=c(4.1,4.1,2.1,0.1), las=1, pty="s")
plot.rf(trout, chr=c(2,29, 10,18, 12,16, 14,20, 27,31),
        alternate.chrid=TRUE, main="")



######################################################################
# Figure 11.5
#
#  Plot of estimated recombination fractions (upper left) and
#  LOD scores for a test of r=1/2 (lower right) for all pairs of
#  markers on linkage groups 13, 19, ``12.16'' and 25, for the
#  trout data.  Red indicates linkage; blue indicates no
#  linkage.
######################################################################
par(mar=c(4.1,4.1,1.1,2.6), las=1, pty="s")
plot.rf(trout, chr=c(13,19,"12.16",25), main="")



######################################################################
# Figure 11.6
#
#  The genetic map in the trout data plotted against the
#  map estimated from the genotype data.  For each linkage group, the
#  line on the left is that map provided with the data, and the line on
#  the right is the estimated map; line segments connect the
#  positions for each marker.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1,cex.axis=0.8)
  newmap <- est.map(trout, error.prob=0.01, 
                    map.function="kosambi", verbose=FALSE)
plot.map(trout, newmap, main="", alternate.chrid=TRUE,
         xlab="Linkage group")
nL <- summary(newmap)["overall","length"]
oL <- summary.map(trout)["overall","length"]



######################################################################
# Figure 11.7
#
#  LOD curves from a genome scan by Haley--Knott regression for
#  the trout data, with MCE groups included as additive
#  covariates.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex.axis=0.8)
plot(out, ylab="LOD score", alternate.chrid=TRUE,
     xlab="Linkage group")



######################################################################
# Figure 11.8
#
#  LOD curves from a genome scan by Haley--Knott regression for
#  the trout data, with MCE groups included as additive
#  covariates.  The curves in blue are as in
#  Fig. ???; those in red were calculated
#  controlling for a QTL on linkage group 8.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex.axis=0.8)
plot(out, out.c8, chr = -8, col=c("blue","red"), 
     ylab="LOD score", alternate.chrid=TRUE,
     xlab="Linkage group")



######################################################################
# Figure 11.9
#
#  LOD scores, for linkage group 8, from a two-dimensional,
#  two-QTL genome scan with the trout data.  lod_i is
#  displayed in the upper left triangle; lod_fv1 is displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_i and lod_fv1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, lower="cond-int", chr=8,
     mar1=c(4.1,9.6,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 11.10
#
#  LOD scores, for selected linkage groups, from a
#  two-dimensional, two-QTL genome scan with the trout data.
#  lod_i is displayed in the upper left triangle; lod_fv1 is
#  displayed in the lower right triangle.  In the color scale on the
#  right, numbers to the left and right correspond to lod_i and
#  lod_fv1, respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, lower="cond-int", chr=c(7,12.16,13),
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)),
     xlab="Linkage group")



######################################################################
# Figure 11.11
#
#  Plot of the tth phenotype against two-locus genotypes
#  at four pairs of putative linked QTL, for the trout data.
#  Points in red are imputed genotypes.
######################################################################
par(mfrow=c(2,2),cex=0.6,cex.axis=0.8)
set.seed(43229849)
par(mar=c(4.6,4.1,1.6,1.1))
plot.pxg(trout,marker=mar7, main=paste("LG 7: ", 
         paste(mar7, collapse=" x ")), col=c("blue","green2"))
par(mar=c(4.6,5.1,1.6,0.1))
plot.pxg(trout,marker=mar8, main=paste("LG 8: ", 
         paste(mar8, collapse=" x ")), col=c("blue","green2"))
par(mar=c(4.1,4.1,2.1,1.1))
plot.pxg(trout,marker=mar12.16, main=paste("LG 12.16: ", 
         paste(mar12.16, collapse=" x ")), col=c("blue","green2"))
par(mar=c(4.1,5.1,2.1,0.1))
plot.pxg(trout,marker=mar13, main=paste("LG 13: ", 
         paste(mar13, collapse=" x ")), col=c("blue","green2"))



######################################################################
# Figure 11.12
#
#  LOD curves from a genome scan by Haley--Knott regression for
#  the trout data, with MCE groups included as additive
#  covariates (in black), with MCE groups included as interactive
#  covariates (in blue) and for the test of QTL x MCE
#  interaction (in red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex.axis=0.8)
plot(outi, lod=1:3, ylab="LOD score", alternate.chrid=TRUE,
     xlab="Linkage group")



######################################################################
# Figure 11.13
#
#  LOD curves for selected linkage groups from a genome scan by Haley--Knott regression for
#  the trout data, controlling for two interacting loci on
#  linkage group 8, with MCE groups included as additive covariates (in
#  black), with MCE groups included as interactive covariates (in blue)
#  and for the test of QTL x MCE interaction (in red).
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex.axis=0.8)
plot(outi.aq, lod=1:3, ylab="LOD score", alternate.chrid=TRUE,
     xlab="Linkage group", chr=c(8, 9, "10.18", "12.16", 17, 24))



######################################################################
# Figure 11.14
#
#  Profile LOD curves for a six-QTL model, including two
#  epistatic QTL on linkage group 8 and QTL x MCE interactions
#  for the QTL on linkage groups 9 and 10.18, for the trout
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1), las=1, cex=0.8)
plotLodProfile(qtl, col=c("blue","red",rep("black",4)),
               ylab="Profile LOD score", xlab="Linkage group")



######################################################################
# Figure 11.15
#
#  Two-dimensional profile LOD surface for the pair of
#  interacting QTL on linkage group 8, in the context of a six-QTL
#  model, including QTL x MCE interactions for the QTL on
#  linkage groups 9 and 10.18 and additional loci on linkage groups 17
#  and 24, for the trout data.
######################################################################
par(cex.axis=0.9)
plot(out.ap, contour=1.5,
     mar1=c(4.1,9.6,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))


# end of fig11.R
