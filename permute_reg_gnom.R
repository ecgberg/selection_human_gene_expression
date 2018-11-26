#!/usr/bin/Rscript
require(dplyr)
require(tidyr)
require(foreach)

args <- commandArgs(trailingOnly=TRUE)
bin_file <- args[1]
out_obj <- args[2]

ASE_bins <- read.table(bin_file, header=TRUE, stringsAsFactors = FALSE, fill=TRUE)
# ASE_bins$`X.NA.` <- NULL # FOR gnomAD data
ASE_bins$pli_bin <- factor(ASE_bins$pli_bin)
ASE_bins$pli_bin <- factor(ASE_bins$pli_bin,
                           levels=levels(ASE_bins$pli_bin)[c(2,3,1)])
print('spec filtered')
filt_bins <- ASE_bins %>% mutate(gnom_sites=rowSums(.[grep('X[123456789]', names(.))], na.rm=TRUE)) %>%
                    filter(gnom_sites>0 | is.na(X.NA.))
filt_bins <- filt_bins %>% replace_na(list(X1 = 0, X2=0, X3=0, X4=0, X5=0, X6=0, X7=0, X8=0, X9=0))

run_permute <- function(filt_spec, 
                        out_objname) {
  full_spec <- filt_spec %>% mutate(total_sites=rowSums(.[grep('X', names(.))], na.rm=TRUE))
  full_mod <- lm(adj_J2~total_sites, data=full_spec)
  full_mod <- data.frame(summary(full_mod)$coefficients)
  full_mod$bin <- row.names(full_mod)

  
  count_spectra <- filt_spec %>% select(`X1`:`X9`, adj_J2) # FOR gnomAD data
  
  bin_mod <- lm(adj_J2~., data=count_spectra)
  bin_mod <- data.frame(summary(bin_mod)$coefficients)
  bin_mod$bin <- row.names(bin_mod)

  # PERMUTE EVERYTING
  t_full <- foreach(i=1:100, .combine=rbind.data.frame) %do% {
    print(paste('running total rep', i))
    count_spectra$permut_aJ2 <- sample(count_spectra$adj_J2, replace=FALSE, size=length(count_spectra$adj_J2))
    c.J2_byAFbin <- data.frame(summary(lm(permut_aJ2 ~ .-adj_J2, data=count_spectra))$coefficients)
    c.J2_byAFbin$bin <- factor(row.names(c.J2_byAFbin))
    c.J2_byAFbin$i <- i
    return(c.J2_byAFbin)
  }
  
  
  # PERMUTE ACROSS GENES WITHIN AN INDIVIDUAL
  t_ind <- foreach(i=1:100, .combine=rbind.data.frame) %do% {
    print(paste('running across gene rep', i))
    shuffled_spec <- filt_spec %>% group_by(SUBJECT_ID) %>% arrange(by_group=SUBJECT_ID) %>% ungroup()
    aj <- shuffled_spec$adj_J2
    shuffled_spec <- shuffled_spec %>% group_by(SUBJECT_ID) %>% sample_frac(1) %>% ungroup() %>% arrange(group_by=SUBJECT_ID) %>% select(`X1`:`X9`) 
    shuffled_spec$permut_aJ2 <- aj
    
    c.J2_byAFbin <- data.frame(summary(lm(permut_aJ2 ~ ., data=shuffled_spec))$coefficients)
    c.J2_byAFbin$bin <- factor(row.names(c.J2_byAFbin))
    c.J2_byAFbin$i <- i
    return(c.J2_byAFbin)
  }
  
  
  # PERMUTE ACROSS INDIVIDUALS WITHIN A GENE
  t_gene <- foreach(i=1:100, .combine=rbind.data.frame) %do% {
    print(paste('running across ind rep', i))
    shuffled_spec <- filt_spec %>% group_by(ensembl_gene_id) %>% arrange(by_group=ensembl_gene_id) %>% ungroup()
    aj <- shuffled_spec$adj_J2
    shuffled_spec <- shuffled_spec %>% group_by(ensembl_gene_id) %>% sample_frac(1) %>% ungroup() %>% arrange(group_by=ensembl_gene_id) %>% select(`X1`:`X9`) 
    shuffled_spec$permut_aJ2 <- aj
    
    c.J2_byAFbin <- data.frame(summary(lm(permut_aJ2 ~ ., data=shuffled_spec))$coefficients)
    c.J2_byAFbin$bin <- factor(row.names(c.J2_byAFbin))
    c.J2_byAFbin$i <- i
    return(c.J2_byAFbin)
  }
  
  t_full$run <- 'all'
  t_ind$run <- 'ind'
  t_gene$run <- 'gene'
  
  t_all <- rbind.data.frame(t_full, t_ind, t_gene)
  
  save(full_mod, bin_mod, t_all, file=out_objname)
}

run_permute(filt_spec=filt_bins, out_objname=out_obj)
