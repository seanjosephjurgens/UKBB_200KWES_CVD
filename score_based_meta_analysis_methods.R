#!/usr/bin/env Rscript

score_meta <- function(single_variant=F,
                       study_summary_data_list=study_summary_data_list,
                       score_col_vec=score_col_vec,
                       pval_col_vec=pval_col_vec,
                       score.se_col_vec=score.se_col_vec,
                       est_col_vec=est_col_vec,
                       est.se_col_vec=est.se_col_vec,
                       est_type="logistic",
                       calc_raw_meta_fixedeffects_odds_ratio=F,
                       recalc_variance_vec=recalc_variance_vec,
                       mincarriers_col_vec=mincarriers_col_vec,
                       mincarriers_num_vec=min_carriers_vec,
                       meta_mincarriers_num=10,
                       save_output=TRUE,
                       outfile=outfile,
                       make_figures=TRUE,
                       min_number_of_studies_contributing=1,
                       max_meta_maf=NA
                       ){
library(dplyr)
library(tidyr)
#source('/medpop/afib/sjurgens/Rscripts/association_source.R')

#Reformatting and merging files
cat('\nReformatting and merging files...\n')
meta_data <- NULL
for(i in 1:length(names(study_summary_data_list))){

        summary_data <- study_summary_data_list[[i]]
        name <- names(study_summary_data_list)[i]
        score_col <- score_col_vec[i]
        pval_col <- pval_col_vec[i]
        score.se_col <- score.se_col_vec[i]
        recalc_variance <- recalc_variance_vec[i]
        mincarriers_col <- mincarriers_col_vec[i]
        mincarriers_num <- mincarriers_num_vec[i]

        summary_data <- as.data.frame(summary_data)
        summary_data <- summary_data[order(summary_data[,pval_col]),]
        summary_data <- summary_data[which(summary_data[,mincarriers_col]>=mincarriers_num),]
        if(recalc_variance){
                summary_data$Variance <- (summary_data[,score_col])^2 / qchisq(summary_data[,pval_col],
                                          lower.tail=F, df=1)
        }else{
              	summary_data$Variance <- (summary_data[,score.se_col])^2
        }
	colnames(summary_data) <- paste0(name, "_", colnames(summary_data))
        if(single_variant){
                summary_data$variant.id <- rownames(summary_data)
                if(is.null(meta_data)){
                         meta_data <- summary_data
                }else{
                      	meta_data <- merge(meta_data, summary_data, by="variant.id", all=T)
                }
        }else{
              	summary_data$Gene <- rownames(summary_data)

                if(is.null(meta_data)){
                        meta_data <- summary_data
                }else{
                      	meta_data <- merge(meta_data, summary_data, by="Gene", all=T)
                }
        }
}

# Calculating meta-statistics
cat('\nCalculating meta-statistics using a Score-based meta-analysis approach...\n')
meta_data <- meta_data %>% mutate_all(funs(replace_na(.,0)))
study_names <- names(study_summary_data_list)
meta_data$chr <- 'none'
for(i in 1:length(study_names)){
        class(meta_data[,paste0(study_names[i], "_", 'chr')]) <- 'character'
        meta_data[meta_data[,paste0(study_names[i], "_", 'chr')] =='0', paste0(study_names[i], "_", 'chr')] <- 'none'
        meta_data[meta_data$chr=='none', 'chr'] <- meta_data[meta_data$chr=='none', paste0(study_names[i], "_", 'chr')]
}
if(single_variant){
        meta_data$pos <- 0
        for(i in 1:length(study_names)){
                meta_data[meta_data$pos==0, 'pos'] <- meta_data[meta_data$pos==0, paste0(study_names[i], "_", 'pos')]
        }
}else{
      	meta_data$mpos <- 0
        for(i in 1:length(study_names)){
                meta_data[meta_data$mpos==0, 'mpos'] <- meta_data[meta_data$mpos==0, paste0(study_names[i], "_", 'mpos')]
        }
}
score_meta <- 0
variance_meta <- 0
if(single_variant){
        MAC_meta <- 0
        n.obs_meta <- 0
}else{
      	cMAC_meta <- 0
}
for(i in 1:length(names(study_summary_data_list))){
        score_meta <- score_meta + meta_data[,paste0(study_names[i], "_", score_col_vec[i])]
        variance_meta <- variance_meta + meta_data[,paste0(study_names[i], "_Variance")]
        if(single_variant){
                MAC_meta <- MAC_meta + meta_data[,paste0(study_names[i], "_", mincarriers_col_vec[i])]
                n.obs_meta <- n.obs_meta + meta_data[,paste0(study_names[i], "_n.obs")]
        }else{
              	cMAC_meta <- cMAC_meta + meta_data[,paste0(study_names[i], "_", mincarriers_col_vec[i])]
        }
}
meta_data$Score_Meta <- score_meta
meta_data$Variance_Meta <- variance_meta
if(single_variant){
        meta_data$n.obs_Meta <- n.obs_meta
        meta_data$MAC_Meta <- MAC_meta
        meta_data$Freq_Meta <- meta_data$MAC_Meta / meta_data$n.obs_Meta
        meta_data <- meta_data[meta_data$MAC_Meta>=meta_mincarriers_num,]
        if(!is.na(max_meta_maf)){
                meta_data <- meta_data[meta_data$Freq_Meta<=max_meta_maf,]
        }
}else{
      	meta_data$cMAC_meta <- cMAC_meta
        meta_data <- meta_data[meta_data$cMAC_Meta>=meta_mincarriers_num,]
        colnames(meta_data)[ncol(meta_data)] <- paste0(mincarriers_col_vec[1], "_Meta")
}
meta_data$P_Meta <- pchisq(((meta_data$Score_Meta^2) / meta_data$Variance_Meta), lower.tail=F, df=1)
meta_data <- meta_data[order(meta_data$P_Meta),]

meta_data$N_studies_contributing <- 0
cols <- which(colnames(meta_data) %in% paste0(study_names, "_", mincarriers_col_vec))
meta_data$Direction_in_contributing_studies <- NA
cols2 <- which(colnames(meta_data) %in% paste0(study_names, "_", score_col_vec))
meta_data$Est_in_contributing_studies <- NA
cols3 <- which(colnames(meta_data) %in% paste0(study_names, "_", est_col_vec))
for(i in 1:nrow(meta_data)){
  help <- meta_data[i,cols]
  help[help>1] <- 1
  tot <- sum(help)
  meta_data[i,'N_studies_contributing'] <- tot
  help <- meta_data[i,cols2]
  help2 <- meta_data[i,cols3]
  helphelp <- NULL
  helphelp2 <- NULL
  for(jj in 1:length(help)){
    j <- help[jj]
    k <- help2[jj]
    if(j>0){
      if(!is.null(helphelp)){
        helphelp <- paste0(helphelp, ", ")
        helphelp2 <- paste0(helphelp2, ", ")
      }
      helphelp <- paste0(helphelp, "+")
      if(est_type=="logistic"){
        helphelp2 <- paste0(helphelp2, round(exp(k),2))
      }else{
	helphelp2 <- paste0(helphelp2, round(k,2))
      }
    }else if(j<0){
      if(!is.null(helphelp)){
        helphelp <- paste0(helphelp, ", ")
        helphelp2 <- paste0(helphelp2, ", ")
      }
      helphelp <- paste0(helphelp, "-")
      if(est_type=="logistic"){
        helphelp2 <- paste0(helphelp2, round(exp(k),2))
      }else{
	helphelp2 <- paste0(helphelp2, round(k,2))
      }
    }
  }
  meta_data[i,'Direction_in_contributing_studies'] <- helphelp
  meta_data[i,'Est_in_contributing_studies'] <- helphelp2
}

if(calc_raw_meta_fixedeffects_odds_ratio){
meta_data$Raw_Meta_Estimate <- NA
cat('\nEstimating a raw meta-effect size using fixed-effects meta-analysis...')
library(meta)
meta_data <- meta_data[order(meta_data$P_Meta),]
if(single_variant){
        upto <- nrow(meta_data[meta_data$P_Meta<0.005,])
}else{
      	upto <- nrow(meta_data)
}
for(j in 1:upto){
        carriers_colz <- paste0(study_names, '_', mincarriers_col_vec)
        carriers_numz <- meta_data[j, carriers_colz]
        contributing <- which(carriers_numz>0)

        est_colz <- paste0(study_names, '_', est_col_vec)
        est.se_colz <- paste0(study_names, '_', est.se_col_vec)
        studlabz <- study_names
        ests <- meta_data[j, est_colz]
        est.ses <- meta_data[j, est.se_colz]

        studlabz <- studlabz[contributing]
        ests <- ests[contributing]
        est.ses <- est.ses[contributing]

        estsz <- NULL
        est.sesz <- NULL
        for(i in 1:length(ests)){
                estsz <- c(estsz, ests[1,i])
                est.sesz <- c(est.sesz, est.ses[1,i])
        }
	if(meta_data[j, 'N_studies_contributing']<2){
                meta_data[j, 'Raw_Meta_Estimate'] <- paste0(round(exp(estsz[1]),2), " [", round(exp(estsz[1]-1.96*est.sesz[1]),2), "; ", round(exp(estsz[1]+1.96*est.sesz[1]),2), "]")
        }else{
              	mod <- metagen(TE=estsz, seTE=est.sesz, studlab = studlabz,
                        sm = "MD",
                        level=0.95, level.comb =0.95,
                        comb.fixed = T, comb.random = F,
                        null.effect = 0)
                meta_data[j, 'Raw_Meta_Estimate'] <- paste0(round(exp(mod$TE.fixed),2), " [", round(exp(mod$TE.fixed-1.96*mod$seTE.fixed),2), "; ", round(exp(mod$TE.fixed+1.96*mod$seTE.fixed),2), "]")
        }
}
}

meta_data <- meta_data[meta_data$N_studies_contributing>=min_number_of_studies_contributing,]
n_genes <- nrow(meta_data)
lambda <- (median((qchisq(meta_data$P_Meta, df=1, lower.tail=F)), na.rm=T))/.456
alpha <- 0.05/n_genes
n_significant <- nrow(meta_data[meta_data$P_Meta<alpha,])
if(single_variant){
        cat('\nNumber of variants analyzed in meta analysis is\n')
}else{
      	cat('\nNumber of genes analyzed in meta analysis is\n')
}
cat('\t', n_genes, '\n')
cat('\nTest-wide lambda for this analysis is\n')
cat('\t', round(lambda,2), '\n')
cat('\nBonferroni-corrected alpha for this analysis is\n')
cat('\t', alpha, '\n')
if(single_variant){
        cat('\nNumber of variants reaching significance in this analysis is\n')
}else{
      	cat('\nNumber of genes reaching significance in this analysis is\n')
}
cat('\t', n_significant, ':\n\n')
if(single_variant){
        print(meta_data[1:n_significant, c(1, c((ncol(meta_data)-11):ncol(meta_data)))])
}else{
      	print(meta_data[1:n_significant, c(1, c((ncol(meta_data)-9):ncol(meta_data)))])
}

if(save_output){
        cat('\nSaving results...\n')
        write.table(meta_data, file=outfile, col.names=T, row.names=F, quote=F, sep='\t')
}

if(make_figures){
        cat('\nMaking figures...\n')
        Xin <- F
        Yin <- F
        if('X' %in% meta_data$chr){
                meta_data[meta_data$chr=='X',] <- '23'
                Xin <- T
        }
	if('XY' %in% meta_data$chr){
                meta_data[meta_data$chr=='XY',] <- '23'
                Xin <- T
        }
	if('Y' %in% meta_data$chr){
                meta_data[meta_data$chr=='Y',] <- '24'
                Yin <- T
        }
	chrlabs <- unique(meta_data$chr)
        if(Xin){chrlabs <- chrlabs[-which(chrlabs=="23")]}
        if(Yin){chrlabs <- chrlabs[-which(chrlabs=="24")]}
        class(chrlabs) <- "numeric"
        chrlabs <- chrlabs[order(chrlabs)]
        class(chrlabs) <- "character"
        if(Xin){chrlabs <- c(chrlabs, 'X')}
        if(Yin){chrlabs <- c(chrlabs, 'Y')}
        class(meta_data$chr) <- "numeric"
        most_sig_with_space <- -log10(min(meta_data[,'P_Meta']))+2
        ylimmax <- plyr::round_any(most_sig_with_space, 5, f = ceiling)
        if(ylimmax<10){
                ylimmax <- 10
        }
	manhattan_outfile <- gsub(".tsv", "_manhattan.pdf", outfile)
        library(qqman)
        pdf(manhattan_outfile, width=16)
        par(mar=c(5.1, 6, 4.1, 2.1))
        if(single_variant){
                class(meta_data$pos) <- "numeric"
                manhattan(meta_data, bp="pos", chr="chr", snp="variant.id", p="P_Meta",
                         chrlabs=chrlabs, main="", ylim=c(0, ylimmax),
                         col=c("dodgerblue4", "firebrick4"),
                         genomewideline= -log10(0.05/nrow(meta_data)),
                         suggestiveline = -log10(1/nrow(meta_data)),
                         cex.lab=1.5)
        }else{
              	manhattan(meta_data, bp="mpos", chr="chr", snp="Gene", p="P_Meta",
                         chrlabs=chrlabs, main="", ylim=c(0, ylimmax),
                         col=c("dodgerblue4", "firebrick4"),
                         genomewideline= -log10(0.05/nrow(meta_data)),
                         suggestiveline = -log10(1/nrow(meta_data)),
                         cex.lab=1.5)
        }
	par(family='sans', cex=1)
        if(single_variant){
                legend("topleft", legend=paste0('n variants = ', n_genes, '    '), text.col="firebrick4",text.font=3)
        }else{
              	legend("topleft", legend=paste0('n genes = ', n_genes, '    '), text.col="firebrick4",text.font=3)
        }
	dev.off()

        qq_outfile <- gsub(".tsv", "_qqplot.pdf", outfile)
        pdf(qq_outfile)
        par(mar=c(5.1, 6, 4.1, 2.1))
        qq(meta_data$P_Meta,
                col="firebrick4",
                ylim=c(0,ylimmax), xlim=c(0,8))
        par(family='sans', cex=1)
        legend("topleft", legend=paste0('lambda = ', round(lambda,2), '    '), text.col="firebrick4",text.font=2)
        dev.off()
}
cat('Done.\n')
return(meta_data)
}
