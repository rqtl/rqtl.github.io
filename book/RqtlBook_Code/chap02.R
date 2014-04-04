######################################################################
# Code for "A guide to QTL mapping with R/qtl"
# Karl W Broman and Saunak Sen
#
# Chapter 2: Importing and simulating data
#######################################################################

# Get a list of the functions in R/qtl
library(help=qtl)

######################################################################
# 2.1 Importing data
######################################################################

###################################
# 2.1.1 Comma-delimited files
###################################

# Load the library; get a quick view of the arguments 
# for the read.cross function.
library(qtl)
args(read.cross)

# Code to read in the sample data in the "csv" format
mydata <- read.cross("csv", "", "mydata.csv")

# The following are all equivalent
mydata <- read.cross("csv", , "mydata.csv")
mydata <- read.cross("csv", file="mydata.csv")
mydata <- read.cross(format="csv", file="mydata.csv")
mydata <- read.cross(file="mydata.csv", format="csv")
mydata <- read.cross(file="mydata.csv")

# Illustrating directory names on Mac
mydata <- read.cross("csv", "~/Desktop", "mydata.csv")

# Illustrating directory names in Windows
mydata <- read.cross("csv", "c:/My Data", "mydata.csv")

# Illustrating use of na.strings, genotypes, and alleles arguments
mydata <- read.cross("csv", "", "mydata.csv", na.strings="na",
                     genotypes=c("BB","BC","CC"),
                     alleles=c("B","C"))

# For files with ";" as field separators and 
# "," used in place of "." in numbers
mydata <- read.cross("csv", , "mydata.csv", sep=";", dec=",")

# Use of # to specify comments
mydata <- read.cross("csv", , "mydata.csv", comment.char="#")

# Read a data file in the "csvr" format
mydata <- read.cross("csvr", , "mydata_rot.csv")

# Read data files in the "csvs" format
mydata <- read.cross("csvs", genfile="mydata_gen.csv",
                     phefile="mydata_phe.csv")

# Another example of reading data in the "csvs" format
mydata <- read.cross("csvs", "../Data", "mydata_gen.csv",
                     "mydata_phe.csv")

###################################
# 2.1.2 MapMaker/QTL
###################################

# An example of reading data in the "mm" format
mydata <- read.cross("mm", "../Data", "mydata.raw", 
                     "mydata.maps")

###################################
# 2.1.3 QTL Cartographer
###################################

# An example of reading data in the "qtlcart" format
mydata <- read.cross("qtlcart", "../Data", "mydata.cro",
                     "mydata.map")

###################################
# 2.1.4 Map Manager QTX
###################################

# Read data in the "qtx" format
mydata <- read.cross("qtx", "../Data", "mydata.qtx")

# Read data in the "qtx" format and then later estimate the map
mydata <- read.cross("qtx", "", "mydata.qtx", 
                     estimate.map=FALSE)
themap <- est.map(mydata, error.prob=0.001)
mydata <- replace.map(mydata, themap)

######################################################################
# 2.2 Exporting data
######################################################################

# Save selected chromosomes in the listeria data to a file on
# the Mac desktop
data(listeria)
write.cross(listeria, "csv", "~/Desktop/listeria", c(5, 13))

######################################################################
# 2.3 Example data
######################################################################

# List of data sets in R/qtl
data(package="qtl")

######################################################################
# 2.4 Data summaries
######################################################################

# summary of the listeria data
data(listeria)
summary(listeria)

# summary plot of the listeria data
plot(listeria)

# Individual panels from the plot.cross.
plot.missing(listeria)
plot.map(listeria)
plot.pheno(listeria, 1)
plot.pheno(listeria, 2)

# Small summary bits from a cross object
nind(listeria)
nphe(listeria)
totmar(listeria)
nchr(listeria)
nmar(listeria)

######################################################################
# 2.5 Simulating data
######################################################################

###################################
# 2.5.1 Additive models
###################################

# Access the map10 genetic map and plot it.
data(map10)
plot(map10)

# Extract the genetic map from the listeria data
data(listeria)
listmap <- pull.map(listeria)

# Create a genetic map with one autosome
mapA <- sim.map(200, 11, include.x=FALSE, eq.spacing=TRUE)

# A genetic map with random marker spacings
mapB <- sim.map(rep(100, 20), 10)

# A similar map, but without anchoring the telomeres
mapC <- sim.map(rep(100, 20), 10, anchor.tel=FALSE)

# Four autosomes of length 50, 75, 100, 125 cM
L <- c(50, 75, 100, 125)
mapD <- sim.map(L, L/5+1, eq.spacing=TRUE, include.x=FALSE)

# Summary of mapD
summary(mapD)

# simulate a backcross with no QTL
simA <- sim.cross(map10, n.ind=100, type="bc")

# Simulate a backcross with two QTL
a <- 2 * sqrt(0.08 / (1 - 2 * 0.08))
mymodel <- rbind(c(1, 50, a), c(14, 65, a))
simB <- sim.cross(map10, type="bc", n.ind=200, model=mymodel)

# Simulate an intercross with three QTL
mymodel2 <- rbind(c(3, 40, 0.5, 0), c(3, 65, -0.5, 0),
                  c(4, 5, 0.5, 0.5))
simC <- sim.cross(map10, type="f2", n.ind=250, model=mymodel2)

# Simulate a backcross with some missing data and some genotyping errors
simD <- sim.cross(map10, type="bc", n.ind=200, model=mymodel,
                  error.prob=0.01, missing.prob=0.05)

# Simulate intercross data; apply pattern of missing data in listeria
data(listeria)
listmap <- pull.map(listeria)
simE <- sim.cross(listmap, type="f2", n.ind=nind(listeria),
                  model=mymodel2)
for(i in 1:nchr(simE))
  simE$geno[[i]]$data[ is.na(listeria$geno[[i]]$data) ] <- NA

# Simulate an intercross with three QTL, now with crossover interference
simF <- sim.cross(map10, type="f2", n.ind=250, model=mymodel2,
                  m=10)

# Simulate an intercross with three QTL, now via the Stahl model
simG <- sim.cross(map10, type="f2", n.ind=250, model=mymodel2,
                  m=10, p=0.1)

###################################
# 2.5.2 More complex models
###################################

# Simulate a backcross with two epistatic QTL
data(map10)
nullmodel <- rbind(c(4, 25, 0), c(5, 45, 0))
episim <- sim.cross(map10, type="bc", n.ind=200, 
                    model=nullmodel)
qtlg <- episim$qtlgeno
wh <- qtlg[,1]==1 & qtlg[,2]==1
episim$pheno[wh, 1] <- episim$pheno[wh, 1] - 1

# Create a binary version of the phenotype and paste it into the data
binphe <- as.numeric(episim$pheno[,1] > 1)
episim$pheno$affected <- binphe

# Simulate sex; create a 3rd phenotype with a QTL x sex interaction
sex <- sample(0:1, nind(episim), replace=TRUE)
phe3 <- rnorm(nind(episim), 0, 1)
phe3[wh & sex==0] <- phe3[wh & sex==0] - 1.5
phe3[wh & sex==1] <- phe3[wh & sex==1] - 0.5
episim$pheno$pheno3 <- phe3
episim$pheno$sex <- sex

######################################################################
# 2.6 Internal data structure
######################################################################

###################################
# 2.6.1 Experimental cross
###################################

# get access to hyper data; look at its class 
data(hyper)
class(hyper)

# elements of the hyper data set
names(hyper)

# hyper phenotypes
hyper$pheno[1:5,]

# names and class of chromosomes
names(hyper$geno)
sapply(hyper$geno, class)

# components of hyper$geno
names(hyper$geno[[3]])
hyper$geno[[3]]$data[91:94,]
hyper$geno[[3]]$map

# Intermediate results stored within hyper$geno components
names(hyper$geno[[3]])
hyper <- calc.genoprob(hyper, step=10, error.prob=0.01)
names(hyper$geno[[3]])
hyper <- sim.geno(hyper, step=10, n.draws=2, error.prob=0.01)
names(hyper$geno[[3]])
hyper <- calc.errorlod(hyper, error.prob=0.01)
names(hyper$geno[[3]])

# est.rf results added to the cross object
names(hyper)
hyper <- est.rf(hyper)
names(hyper)

# est.rf results in detail
hyper$rf[1:4,1:4]

hyper <- clean(hyper)
names(hyper)
names(hyper$geno[[3]])

###################################
# 2.6.2 Genetic map
###################################

# Get access to map10; look at its class
data(map10)
class(map10)

# Names and class of the components of a genetic map object
names(map10)
sapply(map10, class)

# Individual chromosome
map10[[15]]

# end of chap02.R
