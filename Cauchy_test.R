#! Rscript

############################################################
# Cauchy distribution test, adapted from the STAAR authors!!
############################################################

## Cauchy functions
CCT <- function(pvals, weights=NULL, log10p=FALSE, ignore0s=FALSE, ignore1s=FALSE){
  if(log10p){
    pvals <- 10^(-pvals)
  }
  
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    try(pvals <- pvals[-(which(is.na(pvals))), ], silent=T)
    try(pvals <- pvals[-(which(is.na(pvals)))], silent=T)
  }
  if(sum(is.na(pvals)) > 0){
    stop("ERROR: could not evaluate p-value entries properly, NAs remain in the data and function could not remove these.\n")
  }
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one & !(ignore0s | ignore1s)){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero & ignore0s){
    try(pvals <- pvals[-(which(pvals==0)), ], silent=T)
    try(pvals <- pvals[-(which(pvals==0))], silent=T)
  }else if(is.zero & !ignore0s){
    return(0)
  }
  if(is.one & ignore1s){
    warning("There are p-values that are exactly 1! Ignoring those...")
    try(pvals <- pvals[-(which(pvals==1)), ], silent=T)
    try(pvals <- pvals[-(which(pvals==1))], silent=T)  
  }else if(is.one & !ignore1s){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  if(log10p){
    pval <- -log10(pval)
  }
  return(pval)
}
