library(SPAtest)
library(parallel)
library(future.apply)
library(progressr)

library(data.table)
dat <- fread('META_8strata_both_sexes_EUR_v1_run1.tbl', stringsAsFactors = F, data.table=F)
dat <- dat[which(grepl("canonical", dat$MarkerName)), ]
nrow(dat)
#[1] 1225562

for(study in c("UKB", "CCDG", "Geisinger", "TOPMed", "AoU", "MGB", "FOURIER", "DECLARE")){
  message("Busy with ", study)
  inter <- fread(paste0('../', study, '/', study, '_AF_RVAT_STEP2_v1_both_sexes_EUR_AF.regenie'), stringsAsFactors = F, data.table=F)
  inter$P_signed <- sign(inter$BETA)*10^(-inter$LOG10P) 
  inter <- inter[,c("ID", "BETA", "SE", "P_signed", "cMAC", "N_cases", "N_controls")]
  colnames(inter)[c(2:ncol(inter))] <- paste0(study, "__", colnames(inter)[c(2:ncol(inter))])
  rm <- which(duplicated(inter$ID))
  if(length(rm)>0){inter <- inter[-rm, ]}
  dat <- merge(dat, inter, by.x='MarkerName', by.y='ID', all.x=T)
}


## Define a meta-analysis function
compute_P_SPAgc_fast <- function(df, row_index = NA) {
  
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
                       Cutoff.GC = 1.644854, Cutoff.meta = 1.644854
    )))
  }, error = function(e) {
    # Enhanced error message showing the problematic row index and data
    cat("Error in row:", row_index, "\nProblematic data:", df, "\nError message:", e$message, "\n")
    return(NA)
  })
}

# Prepare the data
study_vec <- c("UKB", "CCDG", "Geisinger", "TOPMed", "AoU", "MGB", "FOURIER", "DECLARE")
dat_filt <- dat[,which(gsub("__.*", "", colnames(dat)) %in% study_vec)]
# requires each study has 4 variables in this order, ordered also by study
dat_filt <- dat_filt[,which(gsub(".*__", "", colnames(dat_filt)) %in% c("P_signed", "cMAC", "N_cases", "N_controls"))]
#dat[c(1:20),'P_SPAgc'] <- pbapply::pbapply(dat_filt[c(1:20), ], 1, function(row) {
#  compute_P_SPAgc_fast(df = row)
#})

### Test on one row
df_test <- as.vector(unlist(as.vector(dat_filt[249878,])))
compute_P_SPAgc_fast(df=df_test)
##
# Use parallel processing to speed up the row-wise operation
n_cores <- detectCores() - 1  # Use one less than the total number of cores to avoid overloading
# Set up the parallel plan with multisession (each worker runs in its own process)
plan(multisession, workers = detectCores() - 1)

# Enable progress handling with progressr
handlers(global = TRUE)
handlers("progress")
# Compute adjusted P-values only for tests with nominal P<0.02
dat$P_SPAgc <- dat$`P-value`
dat_filt2 <- dat_filt[dat$`P-value`<0.05, ]
with_progress({
  p <- progressor(along = c(1:nrow(dat_filt2)))  # Set the correct number of steps
  dat[dat$`P-value`<0.05, 'P_SPAgc'] <- future_sapply(1:nrow(dat_filt2), function(i) {
    p(sprintf("x=%g", i))  # Update the progress bar for each iteration
    compute_P_SPAgc_fast(df = dat_filt2[i, , drop = FALSE], row_index = i)
  })
})

cor(dat$P_SPAgc, 
    dat$`P-value`)

cor(dat[dat$`P-value`<0.05,'P_SPAgc'], 
    dat[dat$`P-value`<0.05,'P-value'])

write.table(dat, file='META_8strata_both_sexes_EUR_v1_run1_P_SPAgc.tsv', col.names=T, row.names=F, quote=F, sep='\t')
system('gzip META_8strata_both_sexes_EUR_v1_run1_P_SPAgc.tsv')

dat <- fread('META_8strata_both_sexes_EUR_v1_run1_P_SPAgc.tsv.gz', stringsAsFactors=F, data.table=F)
