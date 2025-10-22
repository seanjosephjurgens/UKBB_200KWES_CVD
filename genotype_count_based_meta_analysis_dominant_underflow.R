#!/usr/bin/env Rscript

##################################################################################################
# Genotype-count based SPA-adjustment for meta-analysis results, to adjust for case-control imbalance 
#### This version assumes a dominant association model (but is also a good approximation for 
####                                                    very rare variants or collapsing rare variant
####                                                    aggregate tests!
##################################################################################################

##### Expects METAL output from regular beta and se-based meta-analysis, with column 'MarkerName' for the variant
##### Expects study-specific sum stats, including columns "ID", "BETA", "SE", "P_signed", "Ncarriers", "N_CASES", "N_CONTROLS"

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)

meta_metal <- args[1]
adjusted_meta_output <- args[2]
p_cutoff_study <- as.numeric(args[3])
p_cutoff_meta <- as.numeric(args[4])
Cutoff.GC <- qnorm(1 - p_cutoff_study/2)
Cutoff.meta <- qnorm(1 - p_cutoff_meta/2)
study_sumstats_vec <- args[5:(length(args))]

if (!require("SPAtest", character.only = TRUE)) {
  install.packages("SPAtest")
}
if (!require("parallel", character.only = TRUE)) {
  install.packages("parallel")
}
if (!require("future.apply", character.only = TRUE)) {
  install.packages("future.apply")
}
if (!require("progressr", character.only = TRUE)) {
  install.packages("progressr")
}
if (!require("progress", character.only = TRUE)) {
  install.packages("progress")
}
library(SPAtest)
library(parallel)
library(future.apply)
library(progressr)
library(progress)

# Reading in METAL output
dat <- data.table::fread(meta_metal, stringsAsFactors = F, data.table=F)
original_cols <- colnames(dat)

# Reading in and processing study-specific sum stats
for(study_num in c(1:length(study_sumstats_vec))){
  study <- study_sumstats_vec[study_num]
  #message("Busy with ", study)
  inter <- data.table::fread(study, stringsAsFactors = F, data.table=F)
  inter$P <- 2 * pnorm(-abs(inter$BETA / inter$SE))
  if(any(inter$P==0)){inter[inter$P==0, 'P'] <- 1e-300}
  inter$P_signed <- sign(inter$BETA)*inter$P 
  inter <- inter[,c("ID", "BETA", "SE", "P_signed", "Ncarriers", "N_CASES", "N_CONTROLS")]
  colnames(inter)[c(2:ncol(inter))] <- paste0("STUDY", study_num, "__", colnames(inter)[c(2:ncol(inter))])
  rm <- which(duplicated(inter$ID))
  if(length(rm)>0){inter <- inter[-rm, ]}
  dat <- merge(dat, inter, by.x='MarkerName', by.y='ID', all.x=T)
}
#print(head(dat))

## Define a meta-analysis function
compute_P_SPAgc_fast <- function(df, row_index = NA, Cutoff.GC, Cutoff.meta) {
  
  tryCatch({
    
    try(df <- as.vector(unlist(df)), silent=TRUE)
    
    # Check for missing study data
    n_studies <- length(df)/4
    rm <- which(is.na(df[seq(1, length(df), by = 4)]))*4-3
    
    if(length(rm)>0){
      df <- df[-c(outer(rm, 0:3, "+"))]
      n_studies <- length(df)/4
      if(n_studies==1){
        return(abs(df[1]))
      }else if(n_studies==0){
        return(NA)
      }
    }
    # define variables
    p_values <- df[seq(1, length(df), by = 4)]
    count_het <- df[seq(2, length(df), by = 4)]
    count_hom <- rep(0, length(count_het))
    count_case <- df[seq(3, length(df), by = 4)]
    count_control <- df[seq(4, length(df), by = 4)]
    #perform count-based SPA
    return(abs(SPAmeta(pvalue.GC = p_values, 
                       GCmat = cbind(count_het, count_hom), 
                       CCsize.GC = cbind(count_case, count_control), 
                       Cutoff.GC = Cutoff.GC, Cutoff.meta = Cutoff.meta
    )))
  }, error = function(e) {
    # Enhanced error message showing the problematic row index and data
    cat("Error in row:", row_index, "\nProblematic data:", df, "\nError message:", e$message, "\n")
    return(NA)
  })
}

# Prepare the data
study_vec <- paste0("STUDY", c(1:(length(study_sumstats_vec))))
#print(study_vec)
dat_filt <- dat[,which(gsub("__.*", "", colnames(dat)) %in% study_vec)]
#print(head(dat_filt))
# requires each study has 4 variables in this order, ordered also by study
dat_filt <- dat_filt[,which(gsub(".*__", "", colnames(dat_filt)) %in% c("P_signed", "Ncarriers", "N_CASES", "N_CONTROLS"))]
#dat[c(1:20),'P_SPAgc'] <- pbapply::pbapply(dat_filt[c(1:20), ], 1, function(row) {
#  compute_P_SPAgc_fast(df = row)
#})
#print(head(dat_filt))

# Use parallel processing to speed up the row-wise operation
n_cores <- detectCores() - 2  # Use one less than the total number of cores to avoid overloading
# Set up the parallel plan with multisession (each worker runs in its own process)
plan(multisession, workers = detectCores() - 2)

# Enable progress handling with progressr
handlers(global = TRUE)
handlers("progress")

# Compute adjusted P-values only for tests with nominal P<p_cutoff_meta
dat$P_new <- dat$P_old <- as.numeric(dat$`P-value`)
dat_filt2 <- dat_filt[dat$P_old<p_cutoff_meta, ]
#print(head(dat_filt2))
with_progress({
  p <- progressor(along = c(1:nrow(dat_filt2)))  # Set the correct number of steps
  dat[dat$P_old<p_cutoff_meta, 'P_new'] <- future_sapply(1:nrow(dat_filt2), function(i) {
    p(sprintf("x=%g", i))  # Update the progress bar for each iteration
    compute_P_SPAgc_fast(df = dat_filt2[i, , drop = FALSE], Cutoff.GC=Cutoff.GC, Cutoff.meta=Cutoff.meta, row_index = i)
  })
})
dat$P_new <- abs(dat$P_new)

#cor(dat$P_SPAgc, 
#    dat$`P-value`)

#cor(dat[dat$`P-value`<p_cutoff_meta,'P_SPAgc'], 
#    dat[dat$`P-value`<p_cutoff_meta,'P-value'])

# Back-corrected of SE based on BETA and new P-value
dat$StdErr_SPAgc <- dat$StdErr
dat[dat$P_old < p_cutoff_meta, "StdErr_SPAgc"] <-
  abs(dat[dat$P_old < p_cutoff_meta, "Effect"]) /
  qnorm(dat[dat$P_old < p_cutoff_meta, "P_new"] / 2, lower.tail = FALSE)
print(head(dat))

# Fix corrected value, taking into account potential underflow issues
dat$P_SPAgc <- as.character(dat$P_new)
if(any(dat$P_new==0)){
  dat[dat$P_new==0, 'P_SPAgc'] <- dat[dat$P_new==0, 'P-value']
  dat[dat$P_new==0, 'StdErr_SPAgc'] <- dat[dat$P_new==0, 'StdErr']
}

dat <- dat[, (which(colnames(dat)%in%c(original_cols, 'P_SPAgc', 'StdErr_SPAgc')))]
write.table(dat, file=adjusted_meta_output, col.names=T, row.names=F, quote=F, sep='\t')
system(paste0('gzip ', adjusted_meta_output))
