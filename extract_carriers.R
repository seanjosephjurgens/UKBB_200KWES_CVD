#! Rscript

source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")

extract_carriers <- function(groupfile, groupings_to_extract, plinkfile, plinkfile_type="bfile", plink_path, collapse=TRUE, canonical=FALSE, max_maf=0.001){

group <- get(load(groupfile))
group$varid <- paste0(group$chr, ":", group$pos, ":", group$ref, ":", group$alt)

if(canonical){
    group <- group[which(grepl("CANONICAL", group$group_id)), ]
}
group <- group[group$group_id %in% groupings_to_extract, ]
group <- group[,c("varid", "alt", "group_id")]
group$varid <- paste0("chr", group$varid)


final <- NULL
num <- 1
length_nums <- length(unique(group$group_id))
for(grouping in unique(group$group_id)){
    cat('busy with number', num, 'out of', length_nums, '...\n')
    num <- num+1
    write.table(group[group$group_id==grouping, 'varid'], file=paste0('varz_', grouping, '_freq', max_maf, '.tsv'), col.names=F, row.names=F, quote=F)
    write.table(group[group$group_id==grouping,c("varid", "alt")], file=paste0('export-allele_', grouping, '_freq', max_maf, '.tsv'), col.names=F, row.names=F, quote=F)
    try(system(paste0(plink_path, ' ',
                  '--', plinkfile_type, '  ', pfile, '  ',
                  '--extract varz_', grouping, '_freq', max_maf, '.tsv  ',
                  '--make-bed --out bfile_', grouping, '_freq', max_maf
    )))
    try(system(paste0(plink_path, ' ',
                  ' --bfile  bfile_', grouping, '_freq', max_maf,
                  ' --max-maf ', max_maf, ' ',
                  ' --export A --export-allele export-allele_', grouping, "_freq", max_maf, '.tsv',
                  ' --out text_', grouping, '_freq', max_maf
    )))
    if(file.exists(paste0('text_', grouping, '_freq', max_maf, '.raw'))){
        library(data.table)
        library(dplyr)
        raw <- fread(paste0('text_', grouping, '_freq', max_maf, '.raw'), stringsAsFactors=F, data.table=F)
        raw <- raw %>% replace(is.na(.), 0)
    
        if(ncol(raw)==6){
            raw[,paste0(grouping)] <- 0
        }else if(ncol(raw)==7){
            raw[,paste0(grouping)] <- raw[,7]
        }else{
            raw[,paste0(grouping)] <- rowSums(raw[,c(7:(ncol(raw)))])
        }
        if(collapse){
            raw[which(raw[,paste0(grouping)]>1), paste0(grouping)] <-1
        }
        colnames(raw)[(ncol(raw))] <- paste0(grouping, "__freq", max_maf) 
        if(is.null(final)){
            raw <- raw[,c(1:6, (ncol(raw)))]
            final <- raw
        }else{
            raw <- raw[,c(1, (ncol(raw)))]
            final <- merge(final, raw, by="FID", all=T)
        }
    }
    try(system(paste0('rm bfile_', grouping, '_freq', max_maf, '.*')))
    try(system(paste0('rm text_', grouping, '_freq', max_maf, '.*')))
    try(system(paste0('rm varz_', grouping, '_freq', max_maf, '.tsv')))
    try(system(paste0('rm export-allele_', grouping, '_freq', max_maf, '.tsv')))
}

return(final)

}
