#!/usr/bin/env Rscript

########################################################################
# Functions for creating sparse kinship matrix and fitting nullmodels
########################################################################

source('GENESIS_adaptation_source.R')


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
