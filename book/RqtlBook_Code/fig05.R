######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap05.R"
#
# Chapter 5: Non-normal phenotypes
#######################################################################

######################################################################
# Figure 5.1
#
#  LOD scores by nonparametric interval mapping for the
#  listeria data.
######################################################################
par(mar=c(4.1,4.1,0.6,0.1))
plot(out.np, ylab="LOD score", alternate.chrid=TRUE)



######################################################################
# Figure 5.2
#
#  LOD scores by nonparametric interval mapping (in blue) and
#  binary trait mapping (in red) for the listeria data.  The
#  binary trait is defined by survival > 250 hr or not.
######################################################################
par(mar=c(4.1,4.1,0.6,0.1))
plot(out.np, out.bin, col=c("blue", "red"), ylab="LOD score",
     alternate.chrid=TRUE)



######################################################################
# Figure 5.3
#
#  LOD scores from the two-part model, with LOD(pi,mu) in
#  black, LOD(pi) in blue and LOD(mu) in red, for the
#  listeria data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.2p, lodcolumn=1:3, ylab="LOD score", 
     alternate.chrid=TRUE)



######################################################################
# Figure 5.5
#
#  LOD scores from the Cox proportional hazards model, using an
#  approach analogous to Haley--Knott regression, for the listeria 
#  data.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
plot(out.cph, ylab="LOD score", alternate.chrid=TRUE)


# end of fig05.R
