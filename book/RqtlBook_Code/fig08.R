######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap08.R"
#
# Chapter 8: Two-dimensional, two-QTL scans
#######################################################################

######################################################################
# Figure 8.1
#
#  LOD scores, for selected chromosomes, from a two-dimensional,
#  two-QTL genome scan with the hyper data.  lod_i is
#  displayed in the upper left triangle; lod_f is displayed in the
#  lower right triangle.  In the color scale on the right, numbers to
#  the left and right correspond to lod_i and lod_f,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, chr=c(1,4,6,7,15),
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.2
#
#  LOD scores, for selected chromosomes, from a two-dimensional,
#  two-QTL genome scan with the hyper data.  lod_i is
#  displayed in the upper left triangle; lod_fv1 is displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_i and lod_fv1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, chr=c(1,4,6,7,15), lower="cond-int",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.3
#
#  LOD scores, for selected chromosomes, from a two-dimensional,
#  two-QTL genome scan with the hyper data.  lod_a is
#  displayed in the upper left triangle; lod_av1 is displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_a and lod_av1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, chr=c(1,4,6,7,15), upper="add", lower="cond-add",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.4
#
#  LOD scores, for chromosome 1, from a two-dimensional, two-QTL
#  genome scan with the hyper data.  lod_av1 is displayed
#  in the upper left triangle; lod_fv1 is displayed in the
#  lower right triangle.  In the color scale on the right, numbers to
#  the left and right correspond to lod_av1 and lod_fv1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out2, chr=1, upper="cond-add", lower="cond-int",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.5
#
#  LOD scores, for chromosome 3, from a two-dimensional, two-QTL
#  genome scan with the hyper data, with values for pairs of
#  positions that are not separated by a marker replaced by 0.
#  lod_i is displayed in the upper left triangle; lod_av1 is
#  displayed in the lower right triangle.  In the color scale on the
#  right, numbers to the left and right correspond to lod_i and
#  lod_av1, respectively.
######################################################################
par(cex.axis=0.9)
plot(clean(out2), chr=3, lower="cond-add",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.6
#
#  The effect of two putative linked QTL on chromosome 3 on the
#  blood pressure phenotype in the hyper data.  Left: a dot plot
#  of the phenotype as a function of marker genotypes, with black dots
#  corresponding to observed genotypes and red dots corresponding to
#  missing (and so imputed) genotypes.  Right: estimated phenotype
#  averages for each of the four two-locus genotype groups, at the
#  inferred locations of the two putative QTL.
######################################################################
set.seed(1626993)
par(mfrow=c(1,2), mar=c(4.1,4.1,0.1,0.6), cex=0.8)
plot.pxg(hyperc3, marker=mar)
axis(side=1, 0.1, paste(mar[1], ":", sep=""), line=-0.5, xpd=TRUE, tick=FALSE)
axis(side=1, 0.1, paste(mar[2], ":", sep=""), line=0.5, xpd=TRUE, tick=FALSE)
par(mar=c(4.1,4.6,0.1,0.1))
effectplot(hyperc3, mname1="3@37.2", mname2="3@44.7",
           ylim=range(pull.pheno(hyperc3,1)), main="")




######################################################################
# Figure 8.7
#
#  Estimated average blood pressure as a function of genotype at
#  loci on chromosomes 1 and 4 (left panel) or chromosomes 6 and 15
#  (right panel), for the hyper data.
######################################################################
hypersub <- sim.geno(subset(hyper, chr=c(1,4,6,15)), step=2.5,
                     error.prob=0.001, n.draws=256)
par(mfrow=c(1,2), mar=c(4.1,4.1,0.1,0.6), cex=0.8)
effectplot(hypersub, mname1="1@68.3", mname2="4@30",
           ylim=c(95, 110), main="")
par(mar=c(4.1,4.6,0.1,0.1))
effectplot(hypersub, mname1="6@60", mname2="15@18",
           ylim=c(95, 110), main="")



######################################################################
# Figure 8.8
#
#  LOD scores, for selected chromosomes, from a two-dimensional,
#  two-QTL genome scan with the nf1 data, for those individuals
#  receiving the NPcis mutation from their mother.  lod_i is
#  displayed in the upper left triangle; lod_fv1 is displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_i and lod_fv1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out.frommom, chr=c(7,15,17), lower="cond-int",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.9
#
#  LOD scores, for selected chromosomes, from a two-dimensional,
#  two-QTL genome scan with the nf1 data, for those individuals
#  receiving the NPcis mutation from their father.  lod_i is
#  displayed in the upper left triangle; lod_av1 is displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_i and lod_av1,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out.fromdad, chr=c(9,19), lower="cond-add",
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.10
#
#  Estimated proportions of affected individuals as a function
#  of genotype at two putative QTL for individuals receiving the
#  NPcis mutation from their mother (left panel) or from their
#  father (right panel), for the nf1 data.
######################################################################
par(mfrow=c(1,2), mar=c(4.1,4.1,0.1,0.6), cex=0.8)
effectplot(nf1.fm, mname1="7@45", mname2="17@3",
           ylim=c(0,1), main="", ylab="Proportion affected")
par(mar=c(4.1,4.6,0.1,0.1))
effectplot(nf1.fd, mname1="9@55.5", mname2="19@0",
           ylim=c(0,1), main="", ylab="Proportion affected")



######################################################################
# Figure 8.11
#
#  Regions in a two-dimensional, two-QTL scan, with the X
#  chromosome included, that require separate treatment.
######################################################################
par(mar=c(0.2,7.6,0.2,7.6))
plot(0, 0, type="n", xaxs="i", yaxs="i", 
     xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
polygon(c(-1,1,1,-1),c(-1,1,-1,-1), lwd=2, xpd=TRUE)
segments(c(0.65,0.65), c(-1,0.65), c(0.65,1), c(0.65,0.65), lwd=2, xpd=TRUE)
text(0.225, -0.275, "A:A")
text(0.825, -0.275, "A:X")
text(0.81, 0.8, "X:X", adj=c(0,1))



######################################################################
# Figure 8.12
#
#  LOD scores from a two-dimensional, two-QTL genome scan with
#  the gutlength data.  lod_i is displayed in the upper left
#  triangle; lod_fv1 is displayed in the lower right triangle.  In
#  the color scale on the right, numbers to the left and right
#  correspond to lod_i and lod_fv1, respectively.
######################################################################
par(cex.axis=0.9)
plot(out.gl, lower="cond-int", alternate.chrid=TRUE,
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))



######################################################################
# Figure 8.13
#
#  Differences in LOD scores calculated with and without the use
#  of covariates, from a two-dimensional, two-QTL genome scan with the
#  gutlength data.  Differences in lod_i are displayed in the
#  upper left triangle; and differences in lod_f are displayed in
#  the lower right triangle.  In the color scale on the right, numbers
#  to the left and right correspond to lod_i and lod_f,
#  respectively.
######################################################################
par(cex.axis=0.9)
plot(out.gl.a - out.gl, allow.neg=TRUE, alternate.chrid=TRUE,
     mar1=c(4.1,4.1,0.1,1.6), mar2=c(4.1,2.1,0.1,2.1),
     layout=list(cbind(1,2), cbind(4.5,1)))


# end of fig08.R
