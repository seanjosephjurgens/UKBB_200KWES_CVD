#!/usr/bin/env Rscript

##########################################################################################################################################################################################################################
# Source functions for gene-based analyses.
#
# This document consists of two parts: 
# 1) Adaptations of functions from the GENESIS package, v2.18.0. 
#    Functions were modified to allow Saddle-Point-Approximation for gene-based tests, as well as some fixes/adaptations for analysis of very large datasets.
# 2) Wrappers to perform gene-based analyses.
#    Functions were created to easily run gene-based analyses in large datasets including a number of additional functionalities
#
# We want to give a huge shoutout to the GENESIS developers!
##########################################################################################################################################################################################################################

# All analyses were run in R v4.0.0
# Versions of relevant packages / dependencies are below:

# Main packages:
#  [1] Matrix_1.2-18        GWASTools_1.34.0     Biobase_2.48.0      
#  [4] GENESIS_2.18.0       GenomicRanges_1.40.0 GenomeInfoDb_1.24.0 
#  [7] IRanges_2.22.1       S4Vectors_0.26.1     BiocGenerics_0.34.0 
# [10] data.table_1.12.8    SeqVarTools_1.26.0   SeqArray_1.28.1     
# [13] gdsfmt_1.24.1       

#  Dependency packages:
#  [1] zoo_1.8-8              tidyselect_1.1.0       purrr_0.3.4           
#  [4] DNAcopy_1.62.0         splines_4.0.0          lattice_0.20-41       
#  [7] vctrs_0.3.0            generics_0.0.2         GWASExactHW_1.01      
# [10] mgcv_1.8-31            blob_1.2.1             survival_3.1-12       
# [13] rlang_0.4.6            pillar_1.4.4           glue_1.4.1            
# [16] DBI_1.1.0              bit64_0.9-7            GenomeInfoDbData_1.2.3
# [19] foreach_1.5.0          lifecycle_0.2.0        zlibbioc_1.34.0       
# [22] MatrixModels_0.4-1     Biostrings_2.56.0      codetools_0.2-16      
# [25] memoise_1.1.0          SparseM_1.78           lmtest_0.9-37         
# [28] quantreg_5.55          broom_0.5.6            Rcpp_1.0.4.6          
# [31] backports_1.1.7        quantsmooth_1.54.0     XVector_0.28.0        
# [34] bit_1.1-15.2           digest_0.6.25          dplyr_1.0.0           
# [37] grid_4.0.0             bitops_1.0-6           sandwich_2.5-1        
# [40] magrittr_1.5           RCurl_1.98-1.2         tibble_3.0.1          
# [43] RSQLite_2.2.0          mice_3.9.0             SNPRelate_1.22.0      
# [46] crayon_1.3.4           tidyr_1.1.0            pkgconfig_2.0.3       
# [49] ellipsis_0.3.1         logistf_1.23           iterators_1.0.12      
# [52] R6_2.4.1               nlme_3.1-147           compiler_4.0.0        
       

library(data.table)
library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(GenomicRanges)
library(Matrix)
library(SPAtest)

#######################################################
# 1) Adaptations of GENESIS analysis functions

SPA_pval_Sean <- function(score.result, nullmod, G, pval.thresh = 1){
	if (!requireNamespace("SPAtest")) stop("package 'SPAtest' must be installed to calculate SPA p-values")
	expit <- function(x){exp(x)/(1+exp(x))}

	# index for variants with Score.pval < pval.thresh
	# only bother running SPA on these variants
	idx <- which(score.result$Score.pval <= pval.thresh)

	# Set colnames in G to match the IDx
	G <- as.data.frame(G)
	colnames(G) <- c(1:ncol(G))
	#cat('SPA is being applied and Dim G is', dim(G),'...\n')
	
	# update columns in score.result
	setnames(score.result, "Score.pval", "SPA.pval")
        if (nrow(score.result) > 0) {
            score.result$SPA.converged <- NA
        } else {
            return(cbind(score.result, data.frame(SPA.converged=logical())))
        }

	if(length(idx) > 0 ){

		# calculate the fitted values  from the null model; 
                # expit(eta) = expit(X\beta + Zb)
                if (nullmod$family$mixedmodel) {
                    mu <- as.vector(expit(nullmod$workingY - nullmod$resid.conditional))
                } else {
                    mu <- as.vector(expit(nullmod$workingY - nullmod$resid.marginal))
                }

		# W is the diagonal of a matrix
		W <- mu*(1-mu)
		WX <- W*nullmod$model.matrix
		XWX.inv <- solve(crossprod(nullmod$model.matrix,WX))

		# original code filtered variants with MAC <= 3 here

		# loop through variants
		for(i in idx){
			# extract the genotypes
			g <- G[,i]
			# get the score
			s <- score.result$Score[i]

			### "flip" g and s if the minor allele is the reference allele; don't do this, as this is only for single vars and not Burden
			##if(mean(g) > 1){ 
			##	g <- 2-g 
			##	s <- -s
			##}
			# identify which elements in g are hom ref
			homRefSet <- which(g == 0)

			# compute adjusted genotype values
			# G - X(X'WX)^{-1}(X'WG)
			g <- as.vector(g - tcrossprod(nullmod$model.matrix, crossprod(crossprod(WX, g), XWX.inv)))

			# compute variance ratio (similar to SAIGE estimator)
			GPG <- score.result$Score.SE[i]^2
			GWG <- sum(W*g^2)
			r <- GPG/GWG

			# get inputs to SPAtest function
			g <- sqrt(r)*g
			qtilde <- as.numeric(s + crossprod(g, mu))

			# compute SPA p-value
			if(length(homRefSet)/length(g) < 0.5){
	        	tmp <- SPAtest:::Saddle_Prob(q = qtilde, mu = mu, g = g, alpha=5e-8, output = "P")
	        }else{
	        	tmp <- SPAtest:::Saddle_Prob_fast(	q = qtilde, mu = mu, g = g, alpha = 5e-8, output = "P",
	        										gNA = g[homRefSet], gNB = g[-homRefSet], 
	        										muNA = mu[homRefSet], muNB = mu[-homRefSet])
	        }

	        # add in results
	        score.result$SPA.converged[i] <- tmp$Is.converge
	        if(tmp$Is.converge){
	        	score.result$SPA.pval[i] <- tmp$p.value
	        }
		}
	}

	return(score.result)
}


#### Impute to zero function which is appropriate for very rare variant collapsing tests. Added check for overflow error
zeroImpute_Sean <- function(geno, freqz) {
        nrowz <- nrow(geno)
        ncolz <- ncol(geno)
        # Check for overflow size
        try(dimz <- nrowz*ncolz, silent=T)
	      # If fail for overflow, use workaround
        if(!is.na(dimz)){
		            if(dimz < 1500000000){
                        miss.idx <- which(is.na(geno))
                        miss.var.idx <- ceiling(miss.idx/nrow(geno))
                        geno[miss.idx] <- 0
		            }else{
                        for(jk in c(1:ncol(geno))){
                                #cat('Busy with', jk, 'out of', ncol(geno), '...\n')
                                miss.rowz <- which(is.na(geno[,jk]))
                                if(length(miss.rowz)>0){
                                        geno[miss.rowz,jk] <- 0
                                }
                        }

		            }
	      }else{
                for(jk in c(1:ncol(geno))){
                        #cat('Busy with', jk, 'out of', ncol(geno), '...\n')
                        miss.rowz <- which(is.na(geno[,jk]))
                        if(length(miss.rowz)>0){
                                geno[miss.rowz,jk] <- 0
                        }
                }
	      }
	      geno
}

#### Redefine meanImpute function to work for very large numbers of variants (for many variants, there is a glitch that causes it to error)
meanImpute_Sean <- function(geno, freqz) {
	nrowz <- nrow(geno)
	ncolz <- ncol(geno)
	try(dimz <- nrowz*ncolz, silent=T)
	if(!is.na(dimz)){
		if(dimz < 1500000000){
                	miss.idx <- which(is.na(geno))
                	miss.var.idx <- ceiling(miss.idx/nrow(geno))
                	imputed <- 2*freqz[miss.var.idx]
                	geno[miss.idx] <- 2*freqz[miss.var.idx]
		}else{
			for(jk in c(1:ncol(geno))){
				#cat('Busy with', jk, 'out of', ncol(geno), '...\n')
				miss.rowz <- which(is.na(geno[,jk]))
				if(length(miss.rowz)>0){
					imputed <- 2*freqz[jk]
					geno[miss.rowz,jk] <- imputed
				}	
			}
		}
	}else{
                      	for(jk in c(1:ncol(geno))){
                                #cat('Busy with', jk, 'out of', ncol(geno), '...\n')
                                miss.rowz <- which(is.na(geno[,jk]))
                                if(length(miss.rowz)>0){
                                        imputed <- 2*freqz[jk]
                                        geno[miss.rowz,jk] <- imputed
                                }
                        }
	}
	geno
}

## create the burden score, than calls the appropriate single variant test function. 
## can easily implement GxE interaction with the burden score... later!
testVariantSetBurden_Sean <- function(nullmod, G, weights, burden.test, collapse){
    # multiply G by weights and compute burden
    if(is(G, "Matrix")){
        burden <- rowSums(G %*% Diagonal(x = weights))
    }else{
        burden <- colSums(t(G) * weights)
    }
    
    if(collapse){
        burden[which(burden>0)] <- 1
    }

    # adjust burden for covariates and random effects
    Gtilde <- GENESIS:::calcGtilde(nullmod, burden)
    if(is.null(nullmod$RSS0)){
        nullmod$RSS0 <- as.numeric(crossprod(nullmod$Ytilde))
    }
    
    if (burden.test == "Score") {
        out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0) 
    }
    if (burden.test == "Score.SPA") {
	out <- GENESIS:::.testGenoSingleVarScore(Gtilde, G = burden, resid = nullmod$resid, RSS0 = nullmod$RSS0)
        out <- SPA_pval_Sean(score.result = out, nullmod = nullmod, G = burden, pval.thresh = 0.05)
    }
    # if (burden.test == "Wald"){
    #     out <- GENESIS:::.testGenoSingleVarWald(Gtilde, Ytilde = nullmod$Ytilde,
    #                                   n = length(nullmod$Ytilde), k = ncol(nullmod$model.matrix))
    # }
    return(out)
}

## new function that runs both SKAT and fastSKAT
testVariantSetSKAT_Sean <- function(nullmod, G, weights, neig = 200, ntrace = 500, verbose = FALSE){
    # multiply G by weights 
    if(is(G, "Matrix")){
        G <- G %*% Diagonal(x = weights)        
    }else{
        G <- t(t(G) * weights)
    }

    # scores
    U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
    # SKAT test statistic
    Q <- sum(U^2)

    # adjust G for covariates and random effects
    G <- GENESIS:::calcGtilde(nullmod, G) # P^{1/2}GW

    # compute the p-value
    out <- GENESIS:::.calcPvalVCTest(Q = Q, G = G, neig = neig, ntrace = ntrace, verbose = verbose)

    return(list(Q = Q, pval = out$pval, err = out$err, pval.method = out$pval.method))
}

## function for SMMAT and fastSMMAT
testVariantSetSMMAT_Sean <- function(nullmod, G, weights, neig = 200, ntrace = 500, verbose = FALSE) {
    # multiply G by weights 
    if(is(G, "Matrix")){
        G <- G %*% Diagonal(x = weights)        
    }else{
        G <- t(t(G) * weights)
    }

    # scores
    cat('\t\trunning burden component of SMMAT...\n')
    U <- as.vector(crossprod(G, nullmod$resid)) # WGPY
    U.sum <- sum(U) # 1WGPY

    # adjust G for covariates and random effects
    G <- GENESIS:::calcGtilde(nullmod, G) # P^{1/2}GW

    # compute burden p-value
    G.rowSums <- rowSums(G) # P^{1/2}GW1
    GG1 <- crossprod(G, G.rowSums) # WGPGW1  # O(mn)
    V.sum <- sum(GG1) # 1WGPGW1
    burden.pval <- pchisq(U.sum^2/V.sum, df=1, lower.tail=FALSE)

    # adjust U and G for burden
    cat('\t\trunning burden adjustment of test statistics for SMMAT...\n')
    U <- U - GG1*U.sum/V.sum # WGPY - WGPGW1 * 1WGPY/(1WGPGW1)
    G <- G - tcrossprod(G.rowSums, GG1)/V.sum # O(mn)

    # SMMAT test statistic
    cat('\t\trunning adjusted SKAT component of SMMAT...\n')
    Q <- sum(U^2)

    ### alternative to part above; seems to be slower from testing; this is how presented in SMMAT paper ###
    # V <- crossprod(G) # WGPGW  # O(m^2n)
    # GG1 <- rowSums(V) # WGPGW1
    # # denominator for burden
    # V.sum <- sum(GG1) # 1WGPGW1
    # # burden p-value
    # burden.pval <- pchisq(U.sum^2/V.sum, df=1, lower.tail=FALSE)
    # # adjust for burden
    # U <- U - GG1*U.sum/V.sum
    # V <- V - tcrossprod(GG1)/V.sum  # O(m^2)

    # compute the p-value for the "adjusted SKAT" part
    out <- GENESIS:::.calcPvalVCTest(Q = Q, G = G, neig = neig, ntrace = ntrace, verbose = verbose)
    theta.pval <- out$pval
    err <- out$err

    # Fisher's method to combine p-values
    cat('\t\trunning Fishers combination for SMMAT...\n')
    smmat.pval <- tryCatch(pchisq(-2*log(burden.pval)-2*log(theta.pval), df=4, lower.tail = FALSE), error = function(e) { NA })
    if(is.na(smmat.pval)) {
        err <- 1
        smmat.pval <- burden.pval
    }
    return(list(pval_burden = burden.pval, pval_theta = theta.pval, pval_SMMAT = smmat.pval, err = err, pval_theta.method = out$pval.method))
}


testVariantSet_Sean <- function( nullmod, G, weights, use.weights=F,
                            test = c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO"),
                            burden.test = c("Score","Score.SPA"), collapse = FALSE, 
			    vc.type = "regular weighted",
                            neig = 200, ntrace = 500, 
                            rho = seq(from = 0, to = 1, by = 0.1)){
                           # pval.method = c("davies", "kuonen", "liu"),
                           # return.scores = FALSE, return.scores.cov = FALSE){

    test <- match.arg(test)
    burden.test <- match.arg(burden.test)
    vc.type <- match.arg(vc.type)
    # pval.method <- match.arg(pval.method)

    G <- GENESIS:::.genoAsMatrix(nullmod, G)

    if (test == "Burden") {
	if(collapse){
		burden.type <- "collapsing test"
	}else if(!use.weights){
		burden.type <- "regular burden"
	}else{
		burden.type <- "externally weighted burden"
	}
        cat('Running Burden test type', burden.type, 'using Pval method ', burden.test, '...\n')
	out <- testVariantSetBurden_Sean(nullmod, G, weights, burden.test = burden.test, collapse = collapse)
    }
    if (test == "SKAT") {
	cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
                                   # return.scores, return.scores.cov)
    }
    if(test == "fastSKAT"){
        cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if (test == "SMMAT") {
        cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig = Inf, ntrace = Inf)
    }
    if(test == "fastSMMAT"){
        cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSMMAT_Sean(nullmod, G, weights, neig, ntrace)
    }
    if(test == "SKATO"){
        cat('Running variance component-based test type', test, 'type', vc.type, '...\n')
        out <- testVariantSetSKATO_Sean(nullmod, G, weights, rho)
    }
    return(out)
}

setGeneric("assocTestAggregate_Sean", function(gdsobj, ...) standardGeneric("assocTestAggregate_Sean"))

match.arg_Sean <- function(test) {
    if (length(test) > 1) test <- NULL
    match.arg(test, choices=c("Burden", "SKAT", "fastSKAT", "SMMAT", "fastSMMAT", "SKATO"))
}

setMethod("assocTestAggregate_Sean",
          "SeqVarIterator",
          function(gdsobj, null.model, AF.max=1, MAC.max=Inf, use.weights=F,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT", "fastSKAT", "SMMAT", "SKATO"),
                   burden.test=c("Score", "Score.SPA"), collapse=FALSE, vc.type="regular weighted",
                   # pval.method=c("davies", "kuonen", "liu"),
                   neig = 200, ntrace = 500,
                   rho = seq(from = 0, to = 1, by = 0.1),
                   sparse=TRUE, imputed=FALSE,
                   male.diploid=TRUE, genome.build=c("hg19", "hg38"),
                   verbose=TRUE) {

              # check argument values
              test <- match.arg_Sean(test)
              burden.test <- match.arg(burden.test)
              # pval.method <- match.arg(pval.method)

              # don't use sparse matrices for imputed dosages
              if (imputed) sparse <- FALSE

              # coerce null.model if necessary
              if (sparse) null.model <- GENESIS:::.nullModelAsMatrix(null.model)
              
              # filter samples to match null model
              sample.index <- GENESIS:::.setFilterNullModel(gdsobj, null.model, verbose=verbose)

              # do we need to match on alleles?
              match.alleles <- any(c("ref", "alt") %in% names(mcols(currentRanges(gdsobj))))

              # check ploidy
              if (SeqVarTools:::.ploidy(gdsobj) == 1) male.diploid <- FALSE
              
              # results
              res <- list()
              res.var <- list()
              i <- 1
              n.iter <- length(variantFilter(gdsobj))
              set.messages <- ceiling(n.iter / 100) # max messages = 100
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=match.alleles, expanded=TRUE)
                  
                  if (!imputed) {
                      geno <- expandedAltDosage(gdsobj, use.names=FALSE, sparse=sparse)[sample.index,,drop=FALSE]
                  } else {
                      geno <- imputedDosage(gdsobj, use.names=FALSE)[sample.index,,drop=FALSE]
                  }

                  if (match.alleles) {
                      index <- GENESIS:::.matchAlleles(gdsobj, var.info)
                      var.info <- var.info[index,,drop=FALSE]
                      geno <- geno[,index,drop=FALSE]
                  } else {
                      index <- NULL
                  }

                  # number of non-missing samples
                  # n.obs <- colSums(!is.na(geno))
                  n.obs <- GENESIS:::.countNonMissing(geno, MARGIN = 2)
                  
                  # allele frequency
                  freq <- GENESIS:::.alleleFreq(gdsobj, geno, variant.index=index, sample.index=sample.index,
                                      male.diploid=male.diploid, genome.build=genome.build)
                  
                  # filter monomorphic variants
                  keep <- GENESIS:::.filterMonomorphic(geno, count=n.obs, freq=freq$freq, imputed=imputed)

                  # exclude variants with freq > max & MAC > max
                  keep <-  keep & freq$freq <= AF.max & freq$MAC <= MAC.max
                  if (!all(keep)) {
                      var.info <- var.info[keep,,drop=FALSE]
                      geno <- geno[,keep,drop=FALSE]
                      n.obs <- n.obs[keep]
                      freq <- freq[keep,,drop=FALSE]
                  }

                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- GENESIS:::.weightFromFreq(freq$freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- currentVariants(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      if (!is.null(index)) weight <- weight[index]
                      weight <- weight[keep]
                      
                      weight0 <- is.na(weight) | weight == 0
                      if (any(weight0)) {
                          keep <- !weight0
                          var.info <- var.info[keep,,drop=FALSE]
                          geno <- geno[,keep,drop=FALSE]
                          n.obs <- n.obs[keep]
                          freq <- freq[keep,,drop=FALSE]
                          weight <- weight[keep]
                      }
                  }
                  
                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values, unless it is collapsing test in which case we will impute to zero
                      if(collapse){
		          if (any(n.obs < nrow(geno))) {
			  	      geno <- zeroImpute_Sean(geno, freq$freq)
		      	  }
		      }else{
		          if (any(n.obs < nrow(geno))) {
                          	geno <- meanImpute_Sean(geno, freq$freq)
                          }
		      }

            # do the test
            assoc <- testVariantSet_Sean(null.model, G=geno, use.weights=use.weights, weights=weight, 
                                    test=test, burden.test=burden.test, vc.type=vc.type, collapse=collapse,
                                    neig = neig, ntrace = ntrace,
                                    rho=rho)
                                    # pval.method=pval.method)
            res[[i]] <- cbind(res[[i]], assoc, stringsAsFactors=FALSE)
        }

        if (verbose & n.iter > 1 & i %% set.messages == 0) {
            message(paste("Iteration", i , "of", n.iter, "completed"))
        }
        i <- i + 1
        iterate <- SeqVarTools:::iterateFilter(gdsobj, verbose=F)
    }

    res <- list(results=dplyr::bind_rows(res), variantInfo=res.var)
    GENESIS:::.annotateAssoc(gdsobj, res)
})

# Helper function from TopmedPipeline for making GRanges objects from variant grouping files 	  
aggregateGRangesList <- function(variants) {
    stopifnot(all(c("group_id", "chr", "pos") %in% names(variants)))
    groups <- unique(variants$group_id)
    cols <- setdiff(names(variants), c("group_id", "chr", "pos"))
    GRangesList(lapply(setNames(groups, groups), function(g) {
        x <- variants[variants$group_id == g,]
        gr <- GRanges(seqnames=x$chr, ranges=IRanges(start=x$pos, width=1))
        mcols(gr) <- x[,cols]
        gr
    }))
}

# Other useful code: removing related individuals from dataset
removerelated <- function(one_to_one_to_remove_df, all_ind_col, random = TRUE, keep = NULL, fixed = FALSE){
  
  #' Removes related individuals from dataset.
  #' 
  #' Takes a dataframe with two cols of IDs representing pairs of related individuals from a dataset, as well as a single collumn with IDs from all individuals of that (larger) dataset
  #' Outputs the IDs of unrelated individuals ONLY, throwing out related ones.
  #' 
  #' one_to_one_to_remove_df    =    dataframe of two columns representing related pairs (IDs)
  #' all_ind_col                =    single column (as dataframe) with IDs from the entire dataset (may be larger than just related ones)
  #' random                     =    if 'FALSE' will first remove individuals that are related to multiple others in dataset, before randomly excluding individuals from the remaining pairs (WARNING this option is a lot slower, but will maximize eventual sample size)
  #'                                  if 'TRUE' will go down all pairs of related and randomly remove individuals from pairs of related, if both are still present in the dataset 
  #' keep                       =    single column (as dataframe) with IDs that you want to preferentially keep. If an instance arises that an individual from this collumn can be removed, the partner will be preferentially removed. If both individuals of a pair are in this column, one of the two will be excluded randomly.
  #' fixed			=    if 'TRUE' will always remove the same entry during the one-to-one remove phase given the same one_to_one_to_remove_df. The removed entry will always be the SECOND individual in the one-to-one remove entry (the first individual will be kept). Default is 'FALSE'.
  #' Returns a dataframe column with the IDs that made it through the removal process: unrelated individuals                                         
 
  cat("\n")
  cat("\n")
  cat("\n")
  
  cat("Preparing for removal ...\n")
  options(stringsAsFactors = F)
  library(plyr)
  one_to_one_to_remove_df <- as.data.frame(one_to_one_to_remove_df)
  one_to_one_to_remove_df <- cbind(one_to_one_to_remove_df, 
                                   c(1:nrow(one_to_one_to_remove_df)))
  colnames(one_to_one_to_remove_df)[1:3] <- c("ID1", 
                                              "ID2", 
                                              "Match_Number")
  
  all_ind_col <- as.data.frame(all_ind_col)
  all_ind_col <- as.data.frame(cbind(all_ind_col, rep(1, nrow(all_ind_col))))
  colnames(all_ind_col) <- c("ID", "helper") 
  
  cat(paste0("Starting with ", nrow(all_ind_col), " individuals."))
  
  if(!is.null(keep)){
    keep <- as.data.frame(keep)
    keep <- as.data.frame(cbind(keep, rep(1, nrow(keep))))
    colnames(keep) <- c("ID", "keep")
  }
  
  cat("\n")
  cat("Checking to see whether both individuals of all pairs are present in the data ...\n")   
  # We can remove one-to-one matches from the relatedness df if either of the individuals is not present in our cohort
  ## Merge ind. list with first col of matches (all=F), to keep only one-to-one matches were the first ind of a match is present in the dataset
  colnames(all_ind_col)[1] <- "ID1"
  present_line_1 <- merge(all_ind_col, 
                          one_to_one_to_remove_df, 
                          by="ID1", 
                          all=F)
  
  ## Do the same for the second col of matches, to keep only one-to-one matches were the second ind of a match is present in the dataset
  colnames(all_ind_col)[1] <- "ID2"
  present_line_2 <- merge(all_ind_col, 
                          one_to_one_to_remove_df, 
                          by="ID2", 
                          all=F)
  
  cat("Removing these pairs otherwise ...\n")
  ## Now merge these two df's (all=F), to keep only matches where both of the individuals in a match are present in the dataset
  present_both <- merge(present_line_1, 
                        present_line_2, 
                        by=c("Match_Number", "ID1", "ID2"), 
                        all=F)
  present_both <- present_both[,c(2:3)]
  colnames(all_ind_col)[1] <- "ID"
  
  # Tell user there are no related if there are none
  if(nrow(present_both)==0){
    cat("WARNING: no pairs of related in the dataset.\n")
    all_ind_col <- as.data.frame(all_ind_col[,-2])
    colnames(all_ind_col) <- "ID"
    return(all_ind_col)
  }
  
  # Make a frequency table for how often individuals are found in the relatedness table
  col1 <- present_both$ID1
  col2 <- present_both$ID2
  col_of_relatedness <- c(col1, col2)
  freq_table <- as.data.frame(plyr::count(col_of_relatedness))
  
  if(!is.null(keep)){
    cat("\n")
    cat("Keep == active. Before further removing, relatedness between those specified will be evaluated, as to keep specified individuals where possible.\n")
    cat("\n")
    cat("Removing non-specified individiuals where they are paired with those specified by 'keep' ...\n")
    
    cat("0% ...")    
    for(i in 1:nrow(keep)){
      if((i-1) %% 40 == 0){
        cat(paste0("\b\b\b\b\b\b\b\b", round(((i-1)/nrow(keep))*100,0),"% ..."))
      }
      ones_to_remove <- NULL
      continue_1 <- FALSE
      continue_2 <- FALSE
      if(keep[i,1] %in% present_both[,1]){
        ones_to_remove <- c(ones_to_remove, present_both[which(present_both[,1] == keep[i,1]), 2])
        continue_1 <- TRUE
      }
      if(keep[i,1] %in% present_both[,2]){
        ones_to_remove <- c(ones_to_remove, present_both[which(present_both[,2] == keep[i,1]), 1])
        continue_2 <- TRUE
      }
      
      if(continue_1 | continue_2){
        ones_to_remove <- as.data.frame(cbind(ones_to_remove, rep(0, length(ones_to_remove))))
        colnames(ones_to_remove) <- c("ID", "keep")
        all_ind_col <- merge(all_ind_col, ones_to_remove, by = "ID", all=T)
        all_ind_col[is.na(all_ind_col$keep),"keep"] <- 1
        all_ind_col <- all_ind_col[all_ind_col$keep == 1, ]
        all_ind_col <- all_ind_col[,-3]
        
        
          colnames(ones_to_remove)[1] <- "ID1"
          present_both <- merge(present_both, ones_to_remove, by = "ID1", all.x=T, all.y=F)
          present_both[is.na(present_both$keep),"keep"] <- 1
          present_both <- present_both[present_both$keep == 1, ]
          present_both <- present_both[,-3]
        
        
          colnames(ones_to_remove)[1] <- "ID2"
          present_both <- merge(present_both, ones_to_remove, by = "ID2", all.x=T, all.y=F)
          present_both[is.na(present_both$keep),"keep"] <- 1
          present_both <- present_both[present_both$keep == 1, ]
          present_both <- present_both[,-3]
        
      }
    }
    cat("\b\b\b\b\b\b\b\b100% ...              \n")
    present_both$Match_Num <- c(1:nrow(present_both))
    colnames(all_ind_col)[1] <- "ID1"
    present_line_1 <- merge(all_ind_col, 
                            present_both, 
                            by="ID1", 
                            all=F)
    colnames(all_ind_col)[1] <- "ID2"
    present_line_2 <- merge(all_ind_col, 
                            present_both, 
                            by="ID2", 
                            all=F)
    present_both <- merge(present_line_1, 
                          present_line_2, 
                          by=c("Match_Num", "ID1", "ID2"), 
                          all=F)
    present_both <- present_both[,c(2:3)]
    colnames(all_ind_col)[1] <- "ID"
    
    col1 <- present_both$ID1
    col2 <- present_both$ID2
    col_of_relatedness <- c(col1, col2)
    freq_table <- as.data.frame(plyr::count(col_of_relatedness))
    cat("\n")
    cat("Partners of specified individuals have been discarded. If two specified individuals were in a pair, the second entry was removed.\n")
    cat(paste0(nrow(all_ind_col), " individuals remain.\n"))
  }
  
  if(random == FALSE){
    cat("\n")
    cat("Random == FALSE: Individuals related to multiple other individuals will be removed first.\n")
    freq_table <- freq_table[order(freq_table$freq, decreasing = T), ]
    
    cat("\n")
    cat("Removing individuals related to multiple others ...\n")
    cat("0% ...")
    iteration <- 0
    num_to_do <- nrow(freq_table[freq_table$freq > 1, ])
    while(freq_table[1,2] > 1){
      
      if(iteration %% 20 == 0){
        num_left <- nrow(freq_table[freq_table$freq > 1, ])
        cat(paste0("\b\b\b\b\b\b\b\b", round(((num_to_do - num_left)/num_to_do)*100,0), "% ..."))
      }
      iteration <- iteration + 1
      
      all_ind_col <- all_ind_col[- which(all_ind_col$ID == freq_table[1,1]), ]
      
      if(freq_table[1,1] %in% present_both$ID1){
        present_both <- present_both[- which(present_both$ID1 == freq_table[1,1]), ]
      }
      
      if(freq_table[1,1] %in% present_both$ID2){
        present_both <- present_both[- which(present_both$ID2 == freq_table[1,1]), ]
      }
      
      col1 <- present_both$ID1
      col2 <- present_both$ID2
      col_of_relatedness <- c(col1, col2)
      freq_table <- as.data.frame(plyr::count(col_of_relatedness))
      freq_table <- freq_table[order(freq_table$freq, decreasing = T), ]
      
    }
    cat('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100% ...               \n')
    cat(paste0(nrow(all_ind_col), " individuals left after preferentially excluding those with multiple relatedness partners.\n"))
    cat("\n")
    if(!fixed){
    	cat("Only single occurences left: Now removing randomly between pairs of remaining related ...\n")
    	cat("0% ...")
	for(i in 1:nrow(present_both)){
      		if((i-1) %% 200 == 0){
        		cat(paste0("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", round((i / nrow(present_both))*100,0), "% ..."))
      		}
      		id_to_remove <- present_both[i, sample(2,1)]
      		all_ind_col <- all_ind_col[- which(all_ind_col$ID == id_to_remove), ]
    	}
        cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100% ...            \n")
    }else{
    	cat("Only single occurences left: Now removing the second entry per related pair as fixed=TRUE ...\n")
	cat("0%	...")
	for(i in 1:nrow(present_both)){
		if((i-1) %% 200 == 0){
			cat(paste0("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", round((i / nrow(present_both))*100,0), "% ..."))
		}
		id_to_remove <- present_both[i, 2]
		all_ind_col <- all_ind_col[- which(all_ind_col$ID == id_to_remove), ]
	}
	cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100% ...           \n")
    }
    
    cat("Done.\n")
    cat("\n")
    cat(paste0(nrow(all_ind_col), " individuals remain after all removal steps.\n"))
    all_ind_col <- as.data.frame(all_ind_col[,-2])
    colnames(all_ind_col) <- "ID"
    return(all_ind_col)
  }  
  
  
  if(random == TRUE){
    if(!fixed){
    	cat("\n")
    	cat("Random == TRUE: Per pair of related individuals a random individual will be removed each time until no related pairs remain.\n")
    	cat("\n")
    	cat("Removing randomly between related pairs ...\n")
    	num_to_do <- nrow(present_both)
    	cat("0% ...")
	while(nrow(present_both)>0){
      		num_done <- num_to_do - nrow(present_both)
      		if(num_done %% 200 == 0){
        		cat(paste0("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", round((num_done / num_to_do)*100,0), "% ..."))
      		}
      		remove_random_id <- present_both[1,sample(2,1)]
      
      		all_ind_col <- all_ind_col[- which(all_ind_col$ID == remove_random_id), ]
      
      		if(remove_random_id %in% present_both$ID1){
        		present_both <- present_both[- which(present_both$ID1 == remove_random_id), ]
      		}
      
      		if(remove_random_id %in% present_both$ID2){
        		present_both <- present_both[- which(present_both$ID2 == remove_random_id), ]
      		}
      
    	}
	cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100% ...               \n")
    }else{
	cat("\n")
	cat("Random == TRUE: Per pair of related individuals an individual will be removed each time until no related pairs remain.\n")
	cat("fixed == TRUE: The removed individual will always be the second individual.\n")
	cat("\n")
	cat("Removing the second entry from related pairs ...\n")
        num_to_do <- nrow(present_both)
        cat("0% ...")
	while(nrow(present_both)>0){
                num_done <- num_to_do - nrow(present_both)
                if(num_done %% 200 == 0){
                        cat(paste0("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", round((num_done / num_to_do)*100,0), "% ..."))
                }
                remove_random_id <- present_both[1, 2]

                all_ind_col <- all_ind_col[- which(all_ind_col$ID == remove_random_id), ]

                if(remove_random_id %in% present_both$ID1){
                        present_both <- present_both[- which(present_both$ID1 == remove_random_id), ]
                }

                if(remove_random_id %in% present_both$ID2){
                        present_both <- present_both[- which(present_both$ID2 == remove_random_id), ]
                }

        }
	cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100% ...        \n")
    }

    cat("Done.\n")
    cat("\n")
    cat(paste0(nrow(all_ind_col), " individuals remain after all removal steps.\n"))
    all_ind_col <- as.data.frame(all_ind_col[,-2])
    colnames(all_ind_col) <- "ID"
    return(all_ind_col)
  }
}  

	   


	   
#######################################################
# 2) Wrapper functions to perform gene-based analyses

	   
make_sparse_kinship_matrix <- function(KINGfile, famfile, sparse_cutoff=2^(-9/2), outfile_matrix, compute_unrel=FALSE, relat_cutoff=2^(-9/2), outfile_unrel){
	
	#' KINGfile = string specifying the KING .kin0 file, which is produced when running KING2 to get pair-wise kinship estimates.
	#'	      IDs should match the IDs in the genetic data.
	#' famfile = string specifying the PLINK .fam for the genetic with sample IDs. 
	#' sparse_cutoff = numeric indicating the preferred relatedness cutoff for the sparse matrix. Standard set at 2^(-9/2).
	#' outfile_matrix = string specifying the preferred location for the output kinship matrix file.
	#' compute_unrelated = logical indicating whether to make a .tsv file containing IDs of unrelated samples from the data.
	#' outfile_unrelated = string specifying the preferred location for the output unrelated samples file
	
	king <- fread(KINGfile, stringsAsFactors=F, data.table=F, select=c("ID1", "ID2", "Kinship"))
	colnames(king)[3] <- 'value'
	fam <- fread(famfile, stringsAsFactors=F, data.table=F)

	if(compute_unrel==TRUE){
		relat <- king[king$value >= relat_cutoff,]
		unrel_IDs <- removerelated(relat[,c(2,4)], fam[,1], random=FALSE, fixed=TRUE)
		colnames(unrel_IDs) <- "ID"
		write.table(unrel_IDs, file=outfile_unrel, col.names=T, row.names=F, quote=F, sep='\t')
	}	
	
	# scale by 2 so the diagonal is 1 and all monozygotic pairs are 1
	king$value <- king$value*2
	class(king$ID1) <- 'character'
	class(king$ID2) <- 'character'
	
	pheno <- fam[1:2]
	colnames(pheno)[1] <- "ID1"
	king <- merge(king, pheno, by="ID1", all=F)
	colnames(king)[ncol(king)] <- "ID1new"
	colnames(pheno)[1] <- "ID2"
	king <- merge(king, pheno, by="ID2", all=F)
	colnames(king)[ncol(king)] <- "ID2new"
	king$ID1 <- king$ID1new
	king$ID2 <- king$ID2new
	king <- king[,c(1:3)]
	class(king$ID1) <- 'character'
	class(king$ID2) <- 'character'

	#Make sparse kinship matrix
	sparseMat <- makeSparseMatrix(king, thresh = 2*sparse_cutoff, sample.include = sample.id, diag.value = 1, verbose = TRUE)
	save(sparseMat, file=outfile)
}

fit_nullmodel <- function(phenofile, ID_col, Outcome, IV_Rank_Norm=FALSE, 
			  Fixed_Covars=NULL, Test_Covars=NULL, Test_Covar_P_cutoff=0.05,
			  Model_type=c("gaussian", "binomial"), relfile=NULL, unrelfile=NULL, outfile){
	
	#' phenfile = string specifying the phenotype file; phenotype file should be in .tsv format. 
	#'            Phenotype file should contain sample identifiers (that match those in the genetic data), the outcome variable, and any fixed-effects covariates.
	#' ID_col = string specifying the column name for the column containing the sample ID information
	#' Outcome = string specifying the column name for the column containing the information on the outcome variable
	#' IV_Rank_Norm = logical specifying whether to inverse-rank normalize the outcome variable. Only works for quantitative traits.
	#' Fixed_Covars = vector of strings specifiyig the column names for all the fixed covariates the user wants to definitely include in the model.
	#' Test_Covars = vector of strings specifiyig the column names for all the fixed covariates the user wants to include if beneath a certain significance.
	#' Test_Covar_P_cutoff = numeric for P-value cutoff the user wants to use to include tested covariates in the model. Standard is 0.05.
	#' Model_type = string specifying the type of model to use (gaussian for quantitative trait, binomial for binary trait).
	#' relfile = string specifying the kinship matrix file; this file should be a 'dgCMatrix' object saved inside an .RData file. 	      
	#'	     The dgCMatrix object should be based on sample IDs matching the IDs in the genetic data file and phenotype file.
	#'	     This is an optional functionality; this should be used if the user wants to perform 
	#' unrelfile = string specifying a .tsv file containing IDs of samples determined to be unrelated within your dataset.
	#'	       IDs should match the IDs from the genetic data and the phenotype file.
	#' 	       This is an optional functionality, and can be used if the user does not have a kinship matrix file but does have a list of unrelated samples.
	#' outfile = string specifyig where to save the new .RData file containing the GENESIS nullmodel.
	
	#Phenotype file
	phen1<-fread(phenofile,header=F,data.table=F,sep="\t")
	
	# Inverse-rank normalize if specifief and quantitative analysis
	if(IV_Rank_Norm==TRUE & Model_type=="gaussian"){
		phen1[,Outcome]<-qnorm((rank(phen1[,Outcome],na.last="keep")-0.5)/sum(!is.na(phen1[,Outcome])))
	}
	
	if(is.null(relfile) | is.null(unrelfile)){
		cat('\nWarning: kinship file and/or unrelated file have been specified as 'NULL', therefore independence of samples will be assumed at certain steps. This is okay if your samples are all unrelated.\n')
	}
	
	# select unrelated samples and test testable covariates to determine inclusion
	assopcs <- NULL
	if(!is.null(Test_Covars)){
		if(!is.null(unrelfile)){
			unrel <- fread(unrelfile, stringsAsFactors=F, data.table=F)
			unrelphen<-phen1[which(phen1[,IDcol] %in% unrel$ID),]
		}else{
			unrelphen <- phen1
		}
		if(is.null(Fixed_Covars)){
			# Make formula and find associated covariates for tested variables
			form0<-as.formula(paste(Outcome, " ~ ", paste0(Test_Covars, collapse="+")))
		}else{
			form0<-as.formula(paste(Outcome, " ~ ", paste0(Test_Covars, collapse="+"), "+", paste(Fixed_Covars, collapse="+")))
		}
					  
		if(Model_type == "gaussian"){
			sum0<-summary(lm(form0, data=unrelphen))
			assopcs<-names(which(sum0$coefficients[2:(length(Test_Covars)+1),4]<Test_Covar_P_cutoff))
		}else if(Model_type == "binomial"){
			sum0<-summary(glm(form0, family=binomial, data=unrelphen))
			assopcs<-names(which(sum0$coefficients[2:(length(Test_Covars)+1),4]<Test_Covar_P_cutoff))
		}
	}	
	
	#Final covariates
	covs <- unique(c(Fixed_Covars, assopcs))
	cat('\nIncluded covariates:', covs, '\n')			  
	
	#Kinship matrix and run nullmodel
	if(!is.null(relfile)){
		relmat<-get(load(relfile))
		relmat1<-relmat[as.character(phen1[,ID_col]), as.character(phen1[,ID_col])]

		names(phen1)[which(colnames(phen1)==ID_col)]<-"scanID"
		scanAnnot <- ScanAnnotationDataFrame(phen1)
	
		nullmod <- fitNullModel(scanAnnot, outcome = Outcome, covars = covs, cov.mat = relmat1, family=Model_type)
	}else{
		names(phen1)[which(colnames(phen1)==ID_col)]<-"scanID"
		scanAnnot <- ScanAnnotationDataFrame(phen1)
		
		nullmod <- fitNullModel(scanAnnot, outcome = Outcome, covars = covs, cov.mat = NULL, family=Model_type)
	}
					  
	save(nullmod,file=outfile)

}	  


	   
	   
# Burden and collapsing tests, options for SPA
perform_burden_collapse <-function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
				   burden.test=c("Score", "Score.SPA"), collapse=TRUE,
				   AF.max=0.001, MAC.max=Inf, use.weights=FALSE){
	#' 
	#' gdsfile = string specifying the file name of the genetic dataset; dataset should be in SeqArray GDS format 
	#' groupfile = string specifyinng the file name of the grouping file; the grouping file contains information of variants to be included in the analysis:
	#'	       The grouping file should be a single dataframe called 'group' that is saved within a .RData file
	#'	       The dataframe should contain the following columns in this order: varid, group_id, chr, pos, ref, alt. All other columns are optional.
	#'	       Optionally, a column named 'weight' can be added for weighted burden tests.
	#'	       An example of a grouping dataframe bellow:
	#'
	#'	       	             varid        group_id chr       pos ref alt         func Dscore
	#'	       1 1:100007074:CTG:C ENSG00000283761   1 100007074 CTG   C hclof_noflag     NA
	#'	       2 1:100007074:CTG:C ENSG00000117620   1 100007074 CTG   C hclof_noflag     NA
	#'	       3   1:100007098:T:C ENSG00000283761   1 100007098   T   C     missense     26
	#'	       4   1:100007098:T:C ENSG00000117620   1 100007098   T   C     missense     26
	#'	       5   1:100007109:C:T ENSG00000283761   1 100007109   C   T hclof_noflag     NA
	#'	       6   1:100007109:C:T ENSG00000117620   1 100007109   C   T hclof_noflag     NA
	#'	         Dtools    Weight gnomAD_AFR_AMR_EAS_NFE_SAS_POPMAX
	#'	       1     NA 1.0000000                                 0
	#'	       2     NA 1.0000000                                 0
	#'	       3     28 0.9285714                                 0
	#'	       4     28 0.9285714                                 0
	#'	       5     NA 1.0000000                                 0
	#'	       6     NA 1.0000000                                 0
	#'
	#' phenfile = string specifying the phenotype file; phenotype file should be in .tsv format. 
	#' 	      Phenotype file should contain sample identifiers (that match those in the GDS file), the outcome variable, and any fixed-effects covariates.
	#' ID_col = string specifying the column name for the column containing the sample ID information
	#' nullfile = string specifying the null-model file; this file contains the null-model that can be made using the 'fitNullModel' function from GENESIS or using our fit_nullmodel function.
	#' outfile = string specifying the preferred output location for the gene-based results. 
	#' burden.test = string specifying the type of test to perform: Either regular "Score" test or "Score.SPA" test.
	#' collapse = logical specifying whether to perform a simple collapsing test (TRUE) or a regular burden test (FALSE).
	#' AF.max = numeric specifying the maximum allele frequency for including variants in the analysis. Variants with MAF>AF.max will be removed.
	#' MAC.max = numeric specifying the maximum minor allele count for including variants in the analysis. Variants with MAC>MAC.max will be removed.
	#' use.weights = logical indicating whether to use external weights in the burden test. Only works for collapse = FALSE. A column called 'weight' should be included in the grouping file.

	
	if(collapse==TRUE){
		burden.type <- "collapsing test"
	}else if(use.weights==F){
		burden.type <- "regular burden"
	}else{
		burden.type <- "externally weighted burden test"	
	}
	
	cat(paste0('\n\nBurden test type is ', burden.type, ' and Pval method is ', burden.test, '.\n\n\n'))
	
	# Samples
	phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
	names(phen1)[which(colnames(phen1)==ID_col)]<-"sample.id"
	id_int <- FALSE
	if(class(phen1$sample.id)=='integer'){
		id_int <- TRUE
		class(phen1$sample.id) <- 'character'
	}
	samid0<-phen1$sample.id
	
	# Read the GDS file
	gds <- seqOpen(gdsfile, allow.duplicate=T)
	samples <- seqGetData(gds, "sample.id")
	if(id_int){class(samples)<-"character"}
	missamples<-samples[!samples %in% samid0]
	misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
	colnames(misphen)<-names(phen1)
	misphen$sample.id<-missamples
	combphen<-rbind(phen1,misphen)
	rownames(combphen)<-combphen$sample.id
	combphen2<-combphen[samples,]
	if(id_int){class(combphen2$sample.id) <- 'integer'}
	
	# Construct a SeqVarData object
	seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))
	
	# Filter the gdsfile
	seqSetFilter(seqData, sample.id=samid0)
	
	# Annotation file
	annot<-get(load(groupfile))
	annot <- as.data.frame(annot)
	#class(annot$chr) <- "numeric"
	class(annot$pos) <- "numeric"
	
	# Grouping file; add weights if weights are selected
	weights.found<-FALSE
	if(use.weights){
		if(!"weight" %in% colnames(annot)){
			cat("\nWARNING: no column named 'weight' found in the grouping file; no weights will be applied.\n")
			gr<-aggregateGRangesList(annot)
		}else if(collapse){
			cat("\nWARNING: weights for collapsing tests are not currently supported. No weights will be applied. Please set collapse to FALSE.\n")
			gr<-aggregateGRangesList(annot)
		}else{
			#annot <- annot[,c("group_id", "chr", "pos", "ref", "alt", "weight")]
			cat("\nuse.weights=T and 'weight' column found in grouping file; variant weights will be applied.\n")
			gr<-aggregateGRangesList(annot)
			weights.found<-TRUE
		}
	}else{
		gr<-aggregateGRangesList(annot)
	}
	
	# Create the iterator
	iterator <- SeqVarListIterator(seqData, variantRanges=gr)
	
	# Load null model
	nullmod<-get(load(nullfile))
	
	# Perform assocation test; apply weights if provided
	if(weights.found){
		assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test="Burden", burden.test=burden.test, vc.type=NULL, collapse = collapse, verbose=TRUE, use.weights=T, weight.user="weight")
	}else{
		assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test="Burden", burden.test=burden.test, vc.type=NULL, collapse = collapse, verbose=TRUE, use.weight=F)
	}
	
	# Save results
	save(assoc,file=outfile)
	seqClose(gds)
}


# Kernell based gene-based tests			  
kernell_variance_component <- function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
				       AF.max=0.001, MAC.max, use.weights=FALSE, vc.test=c("SKAT", "SKATO", "SMMAT")){
	#' 
	#' gdsfile = string specifying the file name of the genetic dataset; dataset should be in SeqArray GDS format 
	#' groupfile = string specifyinng the file name of the grouping file; the grouping file contains information of variants to be included in the analysis:
	#'	       The grouping file should be a single dataframe called 'group' that is saved within a .RData file
	#'	       The dataframe should contain the following columns in this order: varid, group_id, chr, pos, ref, alt. All other columns are optional.
	#'	       Optionally, a column named 'weight' can be added for weighted burden tests.
	#'	       An example of a grouping dataframe bellow:
	#'
	#'	       	             varid        group_id chr       pos ref alt         func Dscore
	#'	       1 1:100007074:CTG:C ENSG00000283761   1 100007074 CTG   C hclof_noflag     NA
	#'	       2 1:100007074:CTG:C ENSG00000117620   1 100007074 CTG   C hclof_noflag     NA
	#'	       3   1:100007098:T:C ENSG00000283761   1 100007098   T   C     missense     26
	#'	       4   1:100007098:T:C ENSG00000117620   1 100007098   T   C     missense     26
	#'	       5   1:100007109:C:T ENSG00000283761   1 100007109   C   T hclof_noflag     NA
	#'	       6   1:100007109:C:T ENSG00000117620   1 100007109   C   T hclof_noflag     NA
	#'	         Dtools    Weight gnomAD_AFR_AMR_EAS_NFE_SAS_POPMAX
	#'	       1     NA 1.0000000                                 0
	#'	       2     NA 1.0000000                                 0
	#'	       3     28 0.9285714                                 0
	#'	       4     28 0.9285714                                 0
	#'	       5     NA 1.0000000                                 0
	#'	       6     NA 1.0000000                                 0
	#'
	#' phenfile = string specifying the phenotype file; phenotype file should be in .tsv format. 
	#' 	      Phenotype file should contain sample identifiers (that match those in the GDS file), the outcome variable, and any fixed-effects covariates.
	#' ID_col = string specifying the column name for the column containing the sample ID information
	#' nullfile = string specifying the null-model file; this file contains the null-model that can be made using the 'fitNullModel' function from GENESIS or using our fit_nullmodel function.
	#' outfile = string specifying the preferred output location for the gene-based results. 
	#' AF.max = numeric specifying the maximum allele frequency for including variants in the analysis. Variants with MAF>AF.max will be removed.
	#' MAC.max = numeric specifying the maximum minor allele count for including variants in the analysis. Variants with MAC>MAC.max will be removed.
	#' use.weights = logical indicating whether to use external weights in the burden test. Only works for collapse = FALSE. A column called 'weight' should be included in the grouping file.
	#' vc.test = vector of kernell-based tests to perform. 

	
	if("Burden" %in% vc.test){
		stop("Burden type test is not supported by this function. For burden use 'perform_burden_collapse()'. Stopping run.")
	} 
	
	if(use.weights==FALSE){
	        vc.type <- "regular weighted"
	}else{
	      	vc.type <- "externally weighted"
	}
	
	cat(paste0('\n\nVariance component test type is ', vc.type, ' ', vc.test, '.\n\n\n'))
	
	# Samples
	phen1<-fread(phenfile,header=T,data.table=F,sep="\t")
	names(phen1)[which(colnames(phen1)==ID_col)]<-"sample.id"
	id_int <- FALSE
	if(class(phen1$sample.id)=='integer'){
		id_int <- TRUE
		class(phen1$sample.id) <- 'character'
	}
	samid0<-phen1$sample.id
	
	# Read gds file
	gds <- seqOpen(gdsfile, allow.duplicate=T)
	samples <- seqGetData(gds, "sample.id")
	if(id_int){class(samples)<-"character"}
	missamples<-samples[!samples %in% samid0]
	misphen<-data.frame(matrix(NA,nrow=length(missamples),ncol=ncol(phen1)))
	colnames(misphen)<-names(phen1)
	misphen$sample.id<-missamples
	combphen<-rbind(phen1,misphen)
	rownames(combphen)<-combphen$sample.id
	combphen2<-combphen[samples,]
	if(id_int){class(combphen2$sample.id) <- 'integer'}
	
	# Construct a SeqVarData object
	seqData <- SeqVarData(gds, sampleData=AnnotatedDataFrame(combphen2))
	
	# Filter the gdsfile
	seqSetFilter(seqData, sample.id=samid0)
	
	# Annotation file
	annot<-get(load(groupfile))
	annot <- as.data.frame(annot)
	#class(annot$chr) <- "numeric"
	class(annot$pos) <- "numeric"
	
	# Grouping file; add weights if weights are selected
	weights.found<-FALSE
	if(use.weights){
	        if(!"weight" %in% colnames(annot)){
	                cat("\nWARNING: no column named 'weight' found in the grouping file; no weights will be applied.\n")
	                gr<-aggregateGRangesList(annot)
	        }else{
	              	#annot <- annot[,c("group_id", "chr", "pos", "ref", "alt", "weight")]
	                cat("\nuse.weights=T and 'weight' column found in grouping file; variant weights will be applied.\n")
	                gr<-aggregateGRangesList(annot)
	                weights.found<-TRUE
	        }
	}else{
	      	gr<-aggregateGRangesList(annot)
	}
	
	# Create the iterator
	iterator <- SeqVarListIterator(seqData, variantRanges=gr)
	
	# Load null model
	nullmod<-get(load(nullfile))
	
	# Perfrom assocation test; apply weights if provided
	if(weights.found){
	        assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, test=vc.test, burden.test="Score", vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weights=T, weight.user="weight")
	}else{
	      	assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, test=vc.test, burden.test="Score", vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weight=F)
	}
	
	# Save results
	save(assoc,file=outfile)
	seqClose(gds)
}

# Function for reading in the results for one phenotype, particularly handy for when analyses were split by chromosome.
summarydata <- function(files, chrs, thre_cMAC=0, add_col=TRUE, add_col_name="Phenotype", add_col_value=NULL){
	
	#' files = vector of strings indicating the names of all the files belonging to the analysis of one phenotype
	#' chrs = vector of strings or numerics indicating the chromosome numbers belonging to each of the included files
	#' thre_cMAC = numeric indicating the cutoff for inclusion of gene-based results in the final results dataframe. Results with cMAC<thre_cMAC will be removed.
	#' add_col = logical indicating whether to add additional information to the final result dataframe. For example, this could be the name of the analyzed phenotype. 'add_col_value' needs to be specified.
	#' add_col_name = string specifying what to call the additional column. For example, 'Phenotype'.
	#' add_col_value = string or numeric specifying what to fill into the new column. For example, 'Atrial_fibrillation_or_flutter'.
	
	sumres<-NULL
	for (num in 1:length(files)){
		outfile<-files[num]
		chr0<-chrs[num]
		res1<-get(load(outfile))
		sum0<-res1$results
		sum0$chr<-chr0
		varinfo0<-res1$variantInfo
		cMAC0<-NULL
		mpos0<-NULL
		for (genenum in c(1:length(varinfo0))){
			cMAC1<-sum(varinfo0[[genenum]]$MAC)
			mpos1<-mean(varinfo0[[genenum]]$pos)
			cMAC0<-c(cMAC0,cMAC1)
			mpos0<-c(mpos0,mpos1)
		}
		sum0$cMAC<-cMAC0
		sum0$mpos<-mpos0
		sum1<-subset(sum0,cMAC>=thre_cMAC)
		sumres<-rbind(sumres,sum1)
	}
	if(add_col == TRUE & !is.null(add_col_value)){
		sumres[,add_col_name] <- add_col_value
	}
	return(sumres)
}
