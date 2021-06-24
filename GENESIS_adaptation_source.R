#!/usr/bin/env Rscript

############################################################################################################################################################
# Source functions for gene-based analyses. These functions are adaptations on original functions from the GENESIS package, v2.18.0. 
# Therefore a huge shoutout to the GENESIS developers!
# Functions were modified to allow Saddle-Point-Approximation for gene-based tests, as well as some fixes/adaptations for analysis of very large datasets. 
############################################################################################################################################################

# All analyses were run in R v4.0.0
# Versions of relevant packages / dependencies are below:

# Main packages:
# [1] Matrix_1.2-18       SeqVarTools_1.26.0  SeqArray_1.28.1    
# [4] gdsfmt_1.24.1       GWASTools_1.34.0    Biobase_2.48.0     
# [7] BiocGenerics_0.34.0 GENESIS_2.18.0     

#  Dependency packages:
#  [1] zoo_1.8-8              tidyselect_1.1.0       purrr_0.3.4           
#  [4] DNAcopy_1.62.0         splines_4.0.0          lattice_0.20-41       
#  [7] vctrs_0.3.0            generics_0.0.2         GWASExactHW_1.01      
# [10] stats4_4.0.0           mgcv_1.8-31            blob_1.2.1            
# [13] survival_3.1-12        rlang_0.4.6            pillar_1.4.4          
# [16] glue_1.4.1             DBI_1.1.0              bit64_0.9-7           
# [19] GenomeInfoDbData_1.2.3 foreach_1.5.0          lifecycle_0.2.0       
# [22] zlibbioc_1.34.0        Biostrings_2.56.0      MatrixModels_0.4-1    
# [25] codetools_0.2-16       memoise_1.1.0          IRanges_2.22.1        
# [28] SparseM_1.78           GenomeInfoDb_1.24.0    lmtest_0.9-37         
# [31] quantreg_5.55          broom_0.5.6            Rcpp_1.0.4.6          
# [34] backports_1.1.7        quantsmooth_1.54.0     S4Vectors_0.26.1      
# [37] XVector_0.28.0         bit_1.1-15.2           digest_0.6.25         
# [40] dplyr_1.0.0            GenomicRanges_1.40.0   grid_4.0.0            
# [43] bitops_1.0-6           sandwich_2.5-1         magrittr_1.5          
# [46] RCurl_1.98-1.2         RSQLite_2.2.0          tibble_3.0.1          
# [49] mice_3.9.0             SNPRelate_1.22.0       crayon_1.3.4          
# [52] tidyr_1.1.0            pkgconfig_2.0.3        ellipsis_0.3.1        
# [55] data.table_1.12.8      logistf_1.23           iterators_1.0.12      
# [58] R6_2.4.1               nlme_3.1-147           compiler_4.0.0        


library(GENESIS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)
library(Matrix)

### new SPA function

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

