#######################################################################
# scanone_cph.R
# For "A guide to QTL mapping with R/qtl", Karl W Broman and Saunak Sen
#
# scanone via Cox proportional hazards model, using H-K regression
#######################################################################
scanone.cph <-
function(cross, pheno.col=1)
{
  require(survival)
  
  pheno <- pull.pheno(cross, pheno.col)
  cross <- subset(cross, ind = !is.na(pheno))
  pheno <- pheno[!is.na(pheno)]

  if(class(pheno) != "Surv")
    stop("Need the phenotype to be of class \"Surv\".")
  
  chrtype <- sapply(cross$geno, class)
  if(any(chrtype=="X")) {
    warning("Dropping X chromosome.")
    cross <- subset(cross, chr=(chrtype != "X"))
  }
  chr <- names(cross$geno)

  result <- NULL
  for(i in 1:nchr(cross)) {
    if(!("prob" %in% names(cross$geno[[i]]))) {
      warning("First running calc.genoprob.")
      cross <- calc.genoprob(cross)
    }
    p <- cross$geno[[i]]$prob

    # pull out map; drop last column of probabilities
    map <- attr(p, "map")
    p <- p[,,-dim(p)[3],drop=FALSE]
      
    lod <- apply(p, 2, function(a,b)
                 diff(coxph(b ~ a)$loglik)/log(10), pheno)

    z <- data.frame(chr=chr[i], pos=map, lod=lod)

    # special names for rows
    w <- names(map)
    o <- grep("^loc-*[0-9]+", w)
    if(length(o) > 0) # inter-marker locations cited as "c*.loc*"
      w[o] <- paste("c",names(cross$geno)[i],".",w[o],sep="")
    rownames(z) <- w

    result <- rbind(result, z)
  }

  class(result) <- c("scanone", "data.frame")
  result
}

# end of scanone_cph.R
