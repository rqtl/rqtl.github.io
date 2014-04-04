######################################################################
# Code for the figures in "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Note: The code below may rely on other code in the book, included
#       in the file "chap02.R"
#
# Chapter 2: Importing and simulating data
#######################################################################

######################################################################
# Figure 2.1
#
#  Part of a data file in the "csv" format, as it might
#  be viewed in a spreadsheet.
######################################################################

library(qtl)
library(qtlbook)
data(ch3c)
par(mar=c(0.1,1.6,1.1,0.1), cex=0.8)
nx <- 10
ny <- 16
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
     xlim=c(0,nx), ylim=c(ny,0), xaxs="i", yaxs="i")
abline(h=1:ny, v=1:nx)
u <- par("usr")
text(1:nx-0.5, u[4]+diff(u[3:4])/ny*0.5, LETTERS[1:nx], 
     xpd=TRUE)
text(u[1]-diff(u[1:2])/nx*0.1, 1:ny-0.5, 1:ny, 
     xpd=TRUE, adj=c(1,0.5))

# the following pads the phenotype with 0's
y <- as.character(ch3c$pheno[,1])
nc <- nchar(y)
for(i in which(nc<5))
  y[i] <- paste(y[i], paste(rep("0",5-nc[i]),collapse=""), sep="")
y[3] <- "-"
ch3c$pheno[,1] <- y

text(1:3-0.5, 0.5, c("pheno", "sex", "pgm"))
for(i in 1:3)
  text(i-0.5, 4:ny-0.5, ch3c$pheno[1:(ny-3),i])

set.seed(100215)
g <- pull.geno(ch3c)
ch <- rep(names(ch3c$geno),nmar(ch3c))
po <- as.character(unlist(pull.map(ch3c)))

# pad with ".0"
wh <- grep("\\.",po)
po[-wh] <- paste(po[-wh], "0", sep=".")

text(4:nx-0.5, 0.5, colnames(g)[1:(nx-3)])
text(4:nx-0.5, 1.5, ch[1:(nx-3)])
text(4:nx-0.5, 2.5, po[1:(nx-3)])
g[sample(prod(dim(g)),1000)] <- 0
for(i in 4:ny) 
  text(4:nx-0.5, i-0.5, c("-","A","H","B")[g[i, 1:(nx-3)]+1])




######################################################################
# Figure 2.4
#
#  Part of a data file in the "csvr" format, as it might
#  be viewed in a spreadsheet.
######################################################################

data(ch3c)
par(mar=c(0.1,1.6,1.1,0.1), cex=0.8)
nx <- 10
ny <- 16
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
     xlim=c(0,nx), ylim=c(ny,0), xaxs="i", yaxs="i")
abline(h=1:ny, v=1:nx)
u <- par("usr")
text(1:nx-0.5, u[4]+diff(u[3:4])/ny*0.5, LETTERS[1:nx], 
     xpd=TRUE)
text(u[1]-diff(u[1:2])/nx*0.1, 1:ny-0.5, 1:ny, 
     xpd=TRUE, adj=c(1,0.5))

# the following pads the phenotype with 0's
y <- as.character(ch3c$pheno[,1])
nc <- nchar(y)
for(i in which(nc<5))
  y[i] <- paste(y[i], paste(rep("0",5-nc[i]),collapse=""), sep="")
y[3] <- "-"
ch3c$pheno[,1] <- y

text(0.5, 1:3-0.5, c("pheno", "sex", "pgm"))
for(i in 1:3)
  text(4:nx-0.5, i-0.5, ch3c$pheno[1:(nx-3),i])

set.seed(100215)
g <- pull.geno(ch3c)
ch <- rep(names(ch3c$geno),nmar(ch3c))
po <- as.character(unlist(pull.map(ch3c)))

# pad with ".0"
wh <- grep("\\.",po)
po[-wh] <- paste(po[-wh], "0", sep=".")

text(0.5, 4:ny-0.5, colnames(g)[1:(ny-3)])
text(1.5, 4:ny-0.5, ch[1:(ny-3)])
text(2.8, 4:ny-0.5, po[1:(ny-3)], adj=c(1,0.5))
g[sample(prod(dim(g)),1000)] <- 0
for(i in 4:nx) 
  text(i-0.5, 4:ny-0.5, c("-","A","H","B")[g[i, 1:(ny-3)]+1])




######################################################################
# Figure 2.5
#
#  Part of the genotype and phenotype data files for an example
#  of the "csvs" format, as they might be viewed in a
#  spreadsheet.
######################################################################

data(ch3c)
#par(mfrow=c(1,2), cex=0.8)
layout(cbind(1,2), widths=c(4,5))
par(mar=c(0.1,1.6,2.6,0.8), cex=0.8)
nx <- 5
ny <- 13

# the following pads the phenotype with 0's
y <- as.character(ch3c$pheno[,1])
nc <- nchar(y)
for(i in which(nc<5))
  y[i] <- paste(y[i], paste(rep("0",5-nc[i]),collapse=""), sep="")
y[3] <- "-"
ch3c$pheno[,1] <- y
ch3c$pheno$id <- 1:nind(ch3c)

plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
     xlim=c(0,4), ylim=c(ny,0), xaxs="i", yaxs="i")
abline(h=1:ny, v=1:5)
u <- par("usr")
text(1:4-0.5, u[4]+diff(u[3:4])/ny*0.5, LETTERS[1:4], 
     xpd=TRUE)
text(u[1]-diff(u[1:2])/nx*0.1, 1:ny-0.5, 1:ny, 
     xpd=TRUE, adj=c(1,0.5))
text(1:4-0.5, 0.5, c("pheno", "sex", "pgm", "id"))
for(i in 1:4)
  text(i-0.5, 2:ny-0.5, ch3c$pheno[1:(ny-1),i])
mtext(side=3,"Phenotype data file", line=1.5)

par(mar=c(0.1,2.4,2.6,0.1))
plot(0,0,type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
     xlim=c(0,nx), ylim=c(ny,0), xaxs="i", yaxs="i")
mtext(side=3,"Genotype data file", line=1.5)
abline(h=1:ny, v=1:nx)
u <- par("usr")
text(1:nx-0.5, u[4]+diff(u[3:4])/ny*0.5, LETTERS[1:nx], 
     xpd=TRUE)
text(u[1]-diff(u[1:2])/nx*0.1, 1:ny-0.5, 1:ny, 
     xpd=TRUE, adj=c(1,0.5))

text(1-0.5, 0.5, "id")
text(1-0.5, 4:ny-0.5, ch3c$pheno[1:(ny-3),4])

set.seed(100215)
g <- pull.geno(ch3c)
ch <- rep(names(ch3c$geno),nmar(ch3c))
po <- as.character(unlist(pull.map(ch3c)))

# pad with ".0"
wh <- grep("\\.",po)
po[-wh] <- paste(po[-wh], "0", sep=".")

text(2:nx-0.5, 0.5, colnames(g)[1:(nx-3)])
text(2:nx-0.5, 1.5, ch[1:(nx-3)])
text(2:nx-0.5, 2.5, po[1:(nx-3)])
g[sample(prod(dim(g)),1000)] <- 0
for(i in 4:ny) 
  text(2:nx-0.5, i-0.5, c("-","A","H","B")[g[i, 1:(nx-3)]+1])




######################################################################
# Figure 2.8
#
#  A genetic map, with approximately 10 cM marker spacing,
#  modeled after the mouse genome and contained in the map10
#  data set in R/qtl.
######################################################################
par(mar=c(4.1,4.1,0.1,0.1))
data(map10)
plot(map10, main="")


# end of fig02.R
