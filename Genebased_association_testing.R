#!/usr/bin/env Rscript

#################################################################################################
# Functions for running Gene-based association tests
# (These functions are also included in the source functions commit (GENESIS_adaptation_source.R)
##################################################################################################

source('GENESIS_adaptation_source.R')

# Burden and collapsing tests, options for SPA
perform_burden_collapse <-function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
				   burden.test=c("Score", "Score.SPA"), collapse=TRUE,
				   AF.max=0.001, MAC.max=Inf, use.weights=FALSE, weight.beta=c(1,1)){
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
	#' weight.beta = vector of length 2 with the parameters of the beta-distribution to be used for variant weighting based on the MAF in the dataset. 
	#'		 Default is c(1,1) which is equivalent to the uniform distribution.
	#'		 Not compatible with use.weight=T.

	
	if(collapse==TRUE){
	        burden.type <- "collapsing test"
	        weight.beta <- c(1,1)
	        cat("Note: because analysis type is collapsing test, the c(1,1) beta distribution (uniform distribution) will be used.\n")
	}else if(use.weights==F){
	        burden.type <- "regular burden"
	}else{
	      	burden.type <- "externally weighted burden test"
	        weight.beta <- c(1,1)
	        cat("Note: because weights are pre-specified, the c(1,1) beta distribution (uniform distribution) will be used.\n")
	}
	
	cat(paste0('\n\nBurden test type is ', burden.type, ' and pvalue method is ', burden.test, ' with beta distribution of ', paste0("(", weight.beta[1], ",", weight.beta[2], ")"), '.\n\n\n'))
	
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
	        assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test="Burden", burden.test=burden.test, vc.type=NULL, collapse = collapse, verbose=TRUE, use.weights=T, weight.user="weight", weight.beta=c(1,1))
	}else{
	      	assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, MAC.max=MAC.max, test="Burden", burden.test=burden.test, vc.type=NULL, collapse = collapse, verbose=TRUE, use.weight=F, weight.beta=weight.beta)
	}
	
	# Save results
	save(assoc,file=outfile)
	seqClose(gds)
}


# Kernell based gene-based tests			  
kernell_variance_component <- function(gdsfile, groupfile, phenfile, ID_col, nullfile, outfile,
				       AF.max=0.001, MAC.max=Inf, use.weights=FALSE, vc.test=c("Score", "Score.SPA"), test=c("SKAT", "SKATO", "SMMAT", "SKAT_SAIGEGENEplus", "Extract"), SAIGEGENEplus_collapse_threshold=10,weight.beta=c(1,1))){
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

	
	if("Burden" %in% test){
	        stop("Burden type test is not supported by this function. For burden use 'hclofburden()'. Stopping run.")
	}
	
	if(use.weights==F){
	        vc.type <- "regular weighted"
	}else{
	      	vc.type <- "externally weighted"
	        weight.beta <- c(1,1)
	        cat("Note: because weights are pre-specified, the c(1,1) beta distribution (uniform distribution) will be used.\n")
	}
	
	cat(paste0('\n\nVariance component test type is ', vc.type, ' ', test, ' using pvalue method ', vc.test, ' with beta distribution of ', paste0("(", weight.beta[1], ",", weight.beta[2], ")"), '.\n\n\n'))
	
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
	        assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weights=T, weight.user="weight",
	                                         SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=c(1,1))
	}else{
	      	assoc <- assocTestAggregate_Sean(iterator, nullmod, AF.max=AF.max, test=test, vc.test=vc.test, vc.type=vc.type, collapse = FALSE, verbose=TRUE, use.weight=F,
	                                         SAIGEGENEplus_collapse_threshold=SAIGEGENEplus_collapse_threshold, weight.beta=weight.beta)
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
