#!/usr/bin/env Rscript

######################################################################################
# Functions for creating sparse kinship matrix and fitting nullmodels
# (These functions are also icluded in the source commit 'GENESIS_adaptation_source.R')
######################################################################################

source('GENESIS_adaptation_source.R')


make_sparse_kinship_matrix_fromKING <- function(KINGfile, famfile, sparse_cutoff=2^(-9/2), outfile_matrix, compute_unrel=FALSE, relat_cutoff=2^(-9/2), outfile_unrel){
	
	#' KINGfile = string specifying the KING .kin0 file, which is produced when running KING2 to get pair-wise kinship estimates.
	#'	      IDs should match the IDs in the genetic data.
	#' famfile = string specifying the PLINK .fam for the genetic data with sample IDs. 
	#' sparse_cutoff = numeric indicating the preferred relatedness cutoff for the sparse matrix. Standard set at 2^(-9/2) on the KING kinship scale.
	#' outfile_matrix = string specifying the preferred location for the output kinship matrix file.
	#' compute_unrel = logical indicating whether to make a .tsv file containing IDs of unrelated samples from the data.
	#' relat_cutoff = numeric specifying the preferred relatedness cutoff for the unrelated sample file output. Value on the KING kinship scale.
	#' outfile_unrel = string specifying the preferred location for the output unrelated samples file
	
	king <- fread(KINGfile, stringsAsFactors=F, data.table=F, select=c("ID1", "ID2", "Kinship"))
	colnames(king)[3] <- 'value'
	fam <- fread(famfile, stringsAsFactors=F, data.table=F)
	sample.id <- fam[,2]
	
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
	save(sparseMat, file=outfile_matrix)
}

make_sparse_kinship_matrix_fromOther <- function(GRfile, ID1_col, ID2_col, GR_col, estimate_scale=1, famfile, sparse_cutoff=2*(2^(-9/2)), outfile_matrix, compute_unrel=FALSE, relat_cutoff=2*(2^(-9/2)), outfile_unrel){
	
	#' GRfile = string specifying a text file with pair-wise genetic relationship/kinship estimates (one pair-wise estimate per row).
	#'	    IDs should match the IDs in the genetic data.
	#' ID1_col = string specifying the column name of the first ID column in the GRfile
	#' ID2_col = string specifying the column name of the second ID column in the GRfile
	#' GR_col = string specifying the column name of the genetic relationship/kinship estimates in the GRfile.
	#' estimate_scale = numeric specifying the scale to apply to the kinship estimates before constructing the matrix.
	#' 		    The scale should be applied in such a way so the diagonal of the matrix (e.g. full relationship) is equal to 1. As an example, for KING kinship estimates a full relationship is equal to 0.5, so a scale of 2 should be applied.
	#'		    Default is 1.
	#' famfile = string specifying the PLINK .fam for the genetic data with sample IDs. 
	#' sparse_cutoff = numeric indicating the preferred relatedness cutoff for the sparse matrix. Standard set at 2*(2^(-9/2)) on the scale where a full genetic relationship is equal to 1.
	#' outfile_matrix = string specifying the preferred location for the output kinship matrix file.
	#' compute_unrel = logical indicating whether to make a .tsv file containing IDs of unrelated samples from the data.
	#' relat_cutoff = numeric specifying the preferred genetic relatedness/kinship cutoff for the unrelated sample file output. Value should be on the the scale where a full genetic relationship is equal to 1.
	#' outfile_unrel = string specifying the preferred location for the output unrelated samples file
	
	GR <- fread(GRfile, stringsAsFactors=F, data.table=F, select=c(ID1_col, ID2_col, GR_col))
	colnames(GR) <- c('ID1', 'ID2', 'value')
	fam <- fread(famfile, stringsAsFactors=F, data.table=F)
	sample.id <- fam[,2]
	
	if(compute_unrel){
		relat <- GR[GR$value >= (relat_cutoff * estimate_scale),]
		unrel_IDs <- removerelated(relat[,c(2,4)], fam[,1], random=FALSE, fixed=TRUE)
		colnames(unrel_IDs) <- "ID"
		write.table(unrel_IDs, file=outfile_unrel, col.names=T, row.names=F, quote=F, sep='\t')
	}	
	
	# scale the estimates
	GR$value <- GR$value * estimate_scale
	class(GR$ID1) <- 'character'
	class(GR$ID2) <- 'character'
	
	pheno <- fam[1:2]
	colnames(pheno)[1] <- "ID1"
	GR <- merge(GR, pheno, by="ID1", all=F)
	colnames(GR)[ncol(GR)] <- "ID1new"
	colnames(pheno)[1] <- "ID2"
	GR <- merge(GR, pheno, by="ID2", all=F)
	colnames(GR)[ncol(GR)] <- "ID2new"
	GR$ID1 <- GR$ID1new
	GR$ID2 <- GR$ID2new
	GR <- GR[,c(1:3)]
	class(GR$ID1) <- 'character'
	class(GR$ID2) <- 'character'

	#Make sparse kinship matrix
	sparseMat <- makeSparseMatrix(GR, thresh = sparse_cutoff, sample.include = sample.id, diag.value = 1, verbose = TRUE)
	save(sparseMat, file=outfile_matrix)
}


fit_nullmodel <- function(phenofile, ID_col, Outcome, IV_Rank_Norm=FALSE, 
			  Fixed_Covars=NULL, Test_Covars=NULL, Test_Covar_P_cutoff=0.05,
			  Model_type=c("gaussian", "binomial"), relfile=NULL, separate.residual.variances=NULL, unrelfile=NULL, outfile){
	
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
	#' separate.residual.variances = A character string specifying the name of a categorical variable in the phenotype file to be used to compute separate residual error variances for heterogeneous groups.
	#'				 For example, this could specify groups from different genetic ancestries or studies.
	#'				 Only applicable for mixed-model setting when analyzing gaussian traits. 
	#' unrelfile = string specifying a .tsv file containing IDs of samples determined to be unrelated within your dataset.
	#'	       IDs should match the IDs from the genetic data and the phenotype file.
	#' 	       This is an optional functionality, and can be used if the user does not have a kinship matrix file but does have a list of unrelated samples.
	#' outfile = string specifyig where to save the new .RData file containing the GENESIS nullmodel.
	
	#Phenotype file
	phen1<-fread(phenofile,header=T,data.table=F,sep="\t")
	
	# Inverse-rank normalize if specifief and quantitative analysis
	if(IV_Rank_Norm & Model_type=="gaussian"){
		phen1[,Outcome]<-qnorm((rank(phen1[,Outcome],na.last="keep")-0.5)/sum(!is.na(phen1[,Outcome])))
	}
	
	if(is.null(relfile) | is.null(unrelfile)){
		cat('\nWarning: kinship file and/or unrelated file have been specified as NULL, therefore independence of samples will be assumed at certain steps. This is okay if your samples are all unrelated.\n')
	}
	
	if(!is.null(unrelfile)){
		unrel <- fread(unrelfile, stringsAsFactors=F, data.table=F)
		unrelphen<-phen1[which(phen1[,ID_col] %in% unrel[,1]),]
	}else{
		unrelphen <- phen1
	}
	
	
	# select unrelated samples and test testable covariates to determine inclusion
	assopcs <- NULL
	if(!is.null(Test_Covars)){
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
	
	#Separate residual variances for certain pre-specified groups
	if(!is.null(separate.residual.variances) & Model_type !="gaussian"){
		cat('Groups for residual variance have been specified, however the model type is not gaussian. Separate residual variances will not be applied.\n')
		separate.residual.variances <- NULL
	}
	if(!is.null(separate.residual.variances) & is.null(relfile)){
		cat('Groups for residual variance have been specified, however no matrix has been provided to compute a mixed-effects model. Separate residual variances will not be applied.\n')
		separate.residual.variances <- NULL
	}
	
	#Final covariates
	covs <- unique(c(Fixed_Covars, assopcs))
	cat('\nIncluded covariates:', covs, '\n')			  
	
	nullmod <- NULL 
	#Kinship matrix and run nullmodel
	if(!is.null(relfile)){
		relmat<-get(load(relfile))
		relmat1<-relmat[as.character(phen1[,ID_col]), as.character(phen1[,ID_col])]

		names(phen1)[which(colnames(phen1)==ID_col)]<-"scanID"
		scanAnnot <- ScanAnnotationDataFrame(phen1)
	
		# Try and run mixed model
		try(nullmod <- fitNullModel(scanAnnot, outcome = Outcome, covars = covs, cov.mat = relmat1, family=Model_type, group.var=separate.residual.variances))
	}
	
	# If not run using mixed-model, or if failed using mixed-model: Run using standard fixed-effects regression model
	if(is.null(nullmod) | class(nullmod)!="GENESIS.nullMixedModel" | T %in% nullmod$zeroFLAG){
		cat("WARNING: Failed fitting mixed nullmodel or nonconvergence mixed-model. Resorting to regular linear regression, using unrelated individuals if provided.\n")
        	names(unrelphen)[which(colnames(unrelphen)==ID_col)]<-"scanID"
		scanAnnot <- ScanAnnotationDataFrame(unrelphen)
        	nullmod <- fitNullModel(scanAnnot, outcome = Outcome, covars = covs, cov.mat = NULL, family=Model_type)
	}
				  
	save(nullmod,file=outfile)

}	  
