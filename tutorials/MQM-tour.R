###################################################
### chunk number 1: setseed
###################################################
set.seed(19696527)


###################################################
### chunk number 2: 
###################################################
library(qtl)
data(map10)
simcross <- sim.cross(map10, type="f2", n.ind=100, missing.prob=0.02)


###################################################
### chunk number 3: missingdata eval=FALSE
###################################################
## geno.image(simcross)


###################################################
### chunk number 4: 
###################################################
geno.image(simcross)


###################################################
### chunk number 5: 
###################################################
# displays warning because MQM ignores the X chromosome in an F2
augmentedcross <- mqmaugment(simcross, minprob=1.0)


###################################################
### chunk number 6: augment1 eval=FALSE
###################################################
## geno.image(augmentedcross)


###################################################
### chunk number 7: 
###################################################
geno.image(augmentedcross)


###################################################
### chunk number 8: 
###################################################
augmentedcross <- mqmaugment(simcross, minprob=0.1)


###################################################
### chunk number 9: augment2 eval=FALSE
###################################################
## geno.image(augmentedcross)


###################################################
### chunk number 10: 
###################################################
geno.image(augmentedcross)


###################################################
### chunk number 11: augment3
###################################################
data(multitrait)
msim5 <- simulatemissingdata(multitrait,5)
msim10 <- simulatemissingdata(multitrait,10)
msim80 <- simulatemissingdata(multitrait,80)


###################################################
### chunk number 12: augment4
###################################################
maug5 <- mqmaugment(msim5)
maug10 <- mqmaugment(msim10,minprob=0.25)
maug80 <- mqmaugment(msim80,minprob=0.80)


###################################################
### chunk number 13: augmentMinProb
###################################################
maug10minprob <- mqmaugment(msim10,minprob=0.001,verbose=TRUE)
maug10minprobImpute <- mqmaugment(msim10,minprob=0.001,strategy="impute",verbose=TRUE)
# check how many individuals are expanded:
nind(maug10minprob)
nind(maug10minprobImpute)


###################################################
### chunk number 14: augment5
###################################################
mqm5 <- mqmscan(maug5)
mqm10 <- mqmscan(maug10)
mqm80 <- mqmscan(maug80)


###################################################
### chunk number 15: augment5b
###################################################
msim5 <- calc.genoprob(msim5)
one5 <- scanone(msim5)
msim10 <- calc.genoprob(msim10)
one10 <- scanone(msim10)
msim80 <- calc.genoprob(msim80)
one80 <- scanone(msim80)


###################################################
### chunk number 16: augment6 eval=FALSE
###################################################
## op <- par(mfrow=c(2,2))
## plot(mqm5,mqm10,mqm80,col=c("green","blue","red"),main="MQM missing data")
## legend("topleft",c("MQM 5%","MQM 10%","MQM 80%"),col=c("green","blue","red"),lwd=1)
## plot(one5,mqm5,main="5% missing", col=c("black", "green"))
## legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)
## plot(one10,mqm10,main="10% missing", col=c("black", "blue"))
## legend("topleft",c("scanone","MQM"),col=c("black","blue"),lwd=1)
## plot(one80,mqm80,main="80% missing", col=c("black", "red"))
## legend("topleft",c("scanone","MQM"),col=c("black","red"),lwd=1)


###################################################
### chunk number 17: 
###################################################
op <- par(mfrow=c(2,2))
plot(mqm5,mqm10,mqm80,col=c("green","blue","red"),main="MQM missing data")
legend("topleft",c("MQM 5%","MQM 10%","MQM 80%"),col=c("green","blue","red"),lwd=1)
plot(one5,mqm5,main="5% missing", col=c("black", "green"))
legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)
plot(one10,mqm10,main="10% missing", col=c("black", "blue"))
legend("topleft",c("scanone","MQM"),col=c("black","blue"),lwd=1)
plot(one80,mqm80,main="80% missing", col=c("black", "red"))
legend("topleft",c("scanone","MQM"),col=c("black","red"),lwd=1)


###################################################
### chunk number 18: 
###################################################
data(multitrait)
maug_min1 <- mqmaugment(multitrait,minprob=1.0)
mqm_min1  <- mqmscan(maug_min1)


###################################################
### chunk number 19: 
###################################################
mgenop <- calc.genoprob(multitrait, step=5)
m_one <- scanone(mgenop)


###################################################
### chunk number 20: 
###################################################
maug <- mqmaugment(multitrait)
mqm  <- mqmscan(maug)


###################################################
### chunk number 21: MinprobMulti
###################################################
plot(m_one, mqm_min1, col=c("black", "green"), lty=1:2)
legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 22: 
###################################################
real_markers <- mqmextractmarkers(mqm)


###################################################
### chunk number 23: 
###################################################
max(mqm)
find.marker(maug,chr=5,pos=35)
multitoset <- find.markerindex(maug,"GH.117C")
setcofactors <- mqmsetcofactors(maug,cofactors=multitoset)
mqm_co1 <- mqmscan(maug, setcofactors)


###################################################
### chunk number 24: Cofactor4multi eval=FALSE
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_co1))
## plot(mqm_co1)


###################################################
### chunk number 25: 
###################################################
# plot after adding first cofactor
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_co1))
plot(mqm_co1)


###################################################
### chunk number 26: Cofactor4bMULTI eval=FALSE
###################################################
## plot(m_one, mqm_co1, col=c("black","green"), lty=1:2)
## legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 27: 
###################################################
plot(m_one, mqm_co1, col=c("black","green"), lty=1:2)
legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 28: 
###################################################
# summary(mqm_co1)
multitoset <- c(multitoset,find.markerindex(maug,find.marker(maug,4,10)))
setcofactors <- mqmsetcofactors(maug,cofactors=multitoset)
mqm_co2 <- mqmscan(maug, setcofactors)


###################################################
### chunk number 29: twowaycomparison eval=FALSE
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_co2))
## plot(mqm_co1, mqm_co2, col=c("blue","green"), lty=1:2)
## legend("topleft",c("one cofactor","two cofactors"),col=c("blue","green"),lwd=1)


###################################################
### chunk number 30: 
###################################################
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_co2))
plot(mqm_co1, mqm_co2, col=c("blue","green"), lty=1:2)
legend("topleft",c("one cofactor","two cofactors"),col=c("blue","green"),lwd=1)


###################################################
### chunk number 31: threewaycomparisonmulti eval=FALSE
###################################################
## plot(mqm, mqm_co1, mqm_co2,
##      col=c("green","red","blue"), lty=1:3)
## legend("topleft",c("no cofactors", "one cofactor","two cofactors"),col=c("green","red","blue"),lwd=1)


###################################################
### chunk number 32: 
###################################################
# plot closeup of threeway comparison
plot(mqm, mqm_co1, mqm_co2,
     col=c("green","red","blue"), lty=1:3)
legend("topleft",c("no cofactors", "one cofactor","two cofactors"),col=c("green","red","blue"),lwd=1)


###################################################
### chunk number 33:  eval=FALSE
###################################################
## autocofactors <- mqmautocofactors(maug,50)
## mqm_auto <- mqmscan(maug, autocofactors)
## setcofactors <- mqmsetcofactors(maug,5)
## mqm_backw <- mqmscan(maug, setcofactors)


###################################################
### chunk number 34: 
###################################################
autocofactors <- mqmautocofactors(maug,50)
mqm_auto <- mqmscan(maug, autocofactors)
setcofactors <- mqmsetcofactors(maug,5)
mqm_backw <- mqmscan(maug, setcofactors)


###################################################
### chunk number 35: ManualAutoStart eval=FALSE
###################################################
## par(mfrow = c(2,1))
## mqmplot.cofactors(maug,autocofactors, justdots=TRUE)
## mqmplot.cofactors(maug,setcofactors, justdots=TRUE)


###################################################
### chunk number 36: 
###################################################
# plot result of cofactor selection
par(mfrow = c(2,1))
mqmplot.cofactors(maug,autocofactors, justdots=TRUE)
mqmplot.cofactors(maug,setcofactors, justdots=TRUE)


###################################################
### chunk number 37: ManualAuto eval=FALSE
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw))
## plot(mqmgetmodel(mqm_auto))


###################################################
### chunk number 38: 
###################################################
# plot result of cofactor backward elimination
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw))
plot(mqmgetmodel(mqm_auto))


###################################################
### chunk number 39: Backward1multi eval=FALSE
###################################################
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw))
## plot(mqm_backw)


###################################################
### chunk number 40: 
###################################################
# plot result of cofactor backward elimination
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw))
plot(mqm_backw)


###################################################
### chunk number 41: Backward2 eval=FALSE
###################################################
## plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
## legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 42: Backward2multi eval=FALSE
###################################################
## plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
## legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 43: 
###################################################
plot(m_one, mqm_backw, col=c("black","green"), lty=1:2)
legend("topleft",c("scanone","MQM"),col=c("black","green"),lwd=1)


###################################################
### chunk number 44: FigLowAlpha eval=FALSE
###################################################
## mqm_backw_low <- mqmscan(maug, setcofactors, cofactor.significance=0.002)
## par(mfrow = c(2,1))
## plot(mqmgetmodel(mqm_backw_low))
## plot(mqm_backw,mqm_backw_low, col=c("blue","green"), lty=1:2)
## legend("topleft",c("Significance=0.02","Significance=0.002"),col=c("blue","green"),lwd=1)


###################################################
### chunk number 45: 
###################################################
mqm_backw_low <- mqmscan(maug, setcofactors, cofactor.significance=0.002)
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_backw_low))
plot(mqm_backw,mqm_backw_low, col=c("blue","green"), lty=1:2)
legend("topleft",c("Significance=0.02","Significance=0.002"),col=c("blue","green"),lwd=1)


###################################################
### chunk number 46: AutoCofactor eval=FALSE
###################################################
## mqmplot.singletrait(mqm_backw_low, extended=TRUE)


###################################################
### chunk number 47: 
###################################################
mqmplot.singletrait(mqm_backw_low, extended=TRUE)


###################################################
### chunk number 48: QTLeffects eval=FALSE
###################################################
## dirresults <- mqmplot.directedqtl(multitrait,mqm_backw_low)


###################################################
### chunk number 49: 
###################################################
dirresults <- mqmplot.directedqtl(multitrait,mqm_backw_low)


###################################################
### chunk number 50: MainEffectsD1 eval=FALSE
###################################################
## plot.pxg(multitrait,marker="GH.117C")


###################################################
### chunk number 51: 
###################################################
plot.pxg(multitrait,marker="GH.117C")


###################################################
### chunk number 52: epistatic1 eval=FALSE
###################################################
## effectplot(multitrait, mname1="GH.117C", mname2="GA1")


###################################################
### chunk number 53: 
###################################################
effectplot(multitrait, mname1="GH.117C", mname2="GA1")


###################################################
### chunk number 54: epistatic2 eval=FALSE
###################################################
## effectplot(multitrait, mname1="PVV4", mname2="GH.117C")


###################################################
### chunk number 55: 
###################################################
effectplot(multitrait, mname1="PVV4", mname2="GH.117C")


###################################################
### chunk number 56: 
###################################################
require(snow)
results <- mqmpermutation(maug,scanfunction=mqmscan,cofactors=setcofactors,n.cluster=2,n.perm=25,batchsize=25)


###################################################
### chunk number 57: 
###################################################
mqmplot.permutations(results)


###################################################
### chunk number 58: 
###################################################

resultsrqtl <- mqmprocesspermutation(results)
summary(resultsrqtl)


###################################################
### chunk number 59: 
###################################################
data(multitrait)
m_imp <- fill.geno(multitrait)
mqmscanfdr(m_imp,mqmscanall,cofactors=setcofactors,n.cluster=2)


###################################################
### chunk number 60: 
###################################################
data(multitrait)
m_imp <- fill.geno(multitrait)
mqm_imp5 <- mqmscan(m_imp,pheno.col=c(1,2,3,4,5),n.cluster=2)


###################################################
### chunk number 61: 
###################################################
mqmplot.multitrait(mqm_imp5,type="image")


###################################################
### chunk number 62: 
###################################################
cofactorlist <- mqmsetcofactors(m_imp,3)
mqm_imp5 <- mqmscan(m_imp,pheno.col=c(1,2,3,4,5),cofactors=cofactorlist,n.cluster=2)


###################################################
### chunk number 63: 
###################################################
mqmplot.multitrait(mqm_imp5,type="image")


###################################################
### chunk number 64: 
###################################################
mqmplot.multitrait(mqm_imp5,type="lines")


###################################################
### chunk number 65: 
###################################################
mqmplot.circle(m_imp, mqm_imp5)


###################################################
### chunk number 66: 
###################################################
mqmplot.circle(m_imp, mqm_imp5, highlight=2)


###################################################
### chunk number 67: 
###################################################
data(locations)
multiloc <- addloctocross(m_imp, locations)
mqmplot.cistrans(mqm_imp5, multiloc, 5, FALSE, TRUE)


###################################################
### chunk number 68: 
###################################################
mqmplot.circle(multiloc,mqm_imp5, highlight=2)


