require(tidyr)
require(foreach)
require(dplyr)

args <- commandArgs(trailingOnly=TRUE)
spec_file <- args[1]
out_obj <- args[2]

print(args)

tissue_spectra <- read.table(spec_file, sep=' ', header=TRUE,
                             stringsAsFactors = FALSE)

tissue_spectra$pli_bin <- factor(tissue_spectra$pli_bin)
tissue_spectra$pli_bin <- factor(tissue_spectra$pli_bin,
                                 levels=levels(tissue_spectra$pli_bin)[c(2,3,1)])
tissue_spectra <- tissue_spectra %>% replace_na(list(X1 = 0, X2=0, X3=0, X4=0, X5=0, X6=0, X7=0, X8=0, X9=0))

print('set spec')
##### HELPER FXNS #####
cfun <- function(list_A, list_B) {
    foreach( A=list_A, B=list_B ) %do% {
      return(rbind.data.frame(A, B)) }
}

# Help jackknife run a little faster
precompute_ids <- function(orig_tab) {
  full_genes <- unique(orig_tab$ensembl_gene_id)
  row_ids <- foreach(gene=full_genes, i=1:length(full_genes)) %do% {
    print(paste(gene, i, 'of', length(full_genes)))
    row_ids=which(orig_tab$ensembl_gene_id==gene)
    return(row_ids)
  }
  names(row_ids) <- full_genes
  return(row_ids)
}

# Run regression 
get_bin_regs <- function(spec) {
  
  count_spectra <- spec %>% select(`X1`:`X9`, adj_J2) 
  print(head(count_spectra))
  
  # run regression, no intercept, order levels
  c.J2_byAFbin <- data.frame(summary(lm(adj_J2 ~ . , data=count_spectra))$coefficients)
  c.J2_byAFbin$bin <- factor(row.names(c.J2_byAFbin))
  
  return(c.J2_byAFbin)
}

# Run a regression to calculate avg effect of each cis-het site on allelic imbalance
get_full_reg <- function(spec) {
  # regress to get avg contribution of each SNP
  spec_counts <-  spec %>% ungroup() %>% mutate(total_sites=rowSums(.[grep('X', names(.))], na.rm=TRUE))
  J2_byTOTAL <-  summary(lm(adj_J2 ~ total_sites, data=spec_counts))
  s.J2_byTOTAL <- data.frame(J2_byTOTAL$coefficients)
  s.J2_byTOTAL$bin=row.names(s.J2_byTOTAL)
  
  return(s.J2_byTOTAL)
}

print('base regresions')

# Get a jacknife sample leaving out the fold-th fold
jack <- function(spec, precomputed_ids, fold) {
  jsamp <- spec[-c(precomputed_ids[[fold]]),]
  return(jsamp)
}

jack_weighted <- function(jack_reg, total_sample, jack_sample, fold) {
      jack_reg$rep <- fold
      jack_reg$fold_size <- total_sample - jack_sample
      jack_reg$fold_ratio <- total_sample / jack_reg$fold_size
      
      # Calculate weighted estimate
      jack_reg$weighted_estimate <- ( 1 - (jack_reg$fold_size / total_sample) ) * jack_reg$Estimate
      
      return(jack_reg)
}

jack_estimates <- function(full_sample,
                           jack_weights,
                           n_jack) {
  # jack_estimate <- ddply(jack_weights, .(bin), summarize, jack_sum=sum(weighted_estimate))
  jack_estimate <- jack_weights %>% group_by(bin) %>% summarize(jack_sum=sum(weighted_estimate))
  jack_estimate <- left_join(full_sample, jack_estimate)
  jack_estimate$jack_est <- ( n_jack * jack_estimate$Estimate ) - jack_estimate$jack_sum
  
  j_vs <- foreach(js=split(jack_weights, jack_weights$bin), .combine=rbind.data.frame) %do% { 
    theta <- jack_estimate[jack_estimate$bin==unique(js$bin),'Estimate']
    
    pseudo_estimate <- theta*js$fold_ratio - ( (js$fold_ratio - 1) * js$Estimate )
    weight_var <- ( pseudo_estimate - jack_estimate[jack_estimate$bin==unique(js$bin),'jack_est'] )**2
    
    js$weighted_variance <- weight_var / (js$fold_ratio - 1)
    return(js)
  }
  
  jack_variance <- j_vs %>% group_by(bin) %>% summarize(jack_var=sum(weighted_variance)/length(unique(rep)))
  jack_estimate <- left_join(jack_estimate, jack_variance)
  
  return(jack_estimate)
}

print('minijacks')

run_jackknife <- function(spec, ids) {
  # Aggregate over the jackknife samples (each with all individual-gene pairs for a single gene left out)
  fullsample_bin_regs <- get_bin_regs(spec)
  fullsample_full_reg <- get_full_reg(spec)
  print('commence')
  
  j <- foreach(i=1:length(unique(spec$ensembl_gene_id)), gene=unique(spec$ensembl_gene_id), .combine=cfun) %do% {
      print(paste0('jack sample: ', gene, ', ', i, ' of ', length(unique(spec$ensembl_gene_id))))
      jack_samp <- jack(spec, precomputed_ids=ids, gene)
      
      jack_full_reg <- get_full_reg(jack_samp)
      wfr <- jack_weighted(jack_reg=jack_full_reg, 
                           total_sample=nrow(spec),
                           jack_sample=nrow(jack_samp),
                           fold=gene)
      
      jack_bin_reg <- get_bin_regs(jack_samp)
      wbr <- jack_weighted(jack_reg=jack_bin_reg, 
                           total_sample=nrow(spec),
                           jack_sample=nrow(jack_samp),
                           fold=gene)
      
      return(list(weight_full_reg=wfr, weight_bin_reg=wbr))
  }
  
  full_jack_estimate <- jack_estimates(full_sample=fullsample_full_reg,
                                       jack_weights = j[[1]],
                                       n_jack = length(unique(spec$ensembl_gene_id)))
  bin_jack_estimate <- jack_estimates(full_sample=fullsample_bin_regs,
                                      jack_weights = j[[2]],
                                      n_jack=length(unique(spec$ensembl_gene_id)))
  
  return(list(full_jack_estimate, bin_jack_estimate))
  
}

print('defined fxns')

# Precompute info
gene_ids <- precompute_ids(orig_tab=tissue_spectra)
# low_pli <- subset(tissue_spectra, 
#                   tissue_spectra$pli_bin=='low')
# high_pli <- subset(tissue_spectra, 
#                    tissue_spectra$pli_bin=='high')
# low_ids <- precompute_ids(orig_tab=low_pli)
# high_ids <- precompute_ids(orig_tab=high_pli)

# Run regressions with jackknife estimate
trial_run <- run_jackknife(spec=tissue_spectra, ids=gene_ids)
full_jack_estimate <- trial_run[[1]] # total het sites vs allelic imbalance
bin_jack_estimate <- trial_run[[2]] # proportion in each bin vs allelic imbalance (after accounting for total hets)

# pli_run <- run_jackknife(spec=low_pli,
#                          ids=low_ids)
# pli_run[[1]]$pli_bin <- 'low'
# pli_run[[2]]$pli_bin <- 'low'
# 
# pli_high <- run_jackknife(spec=high_pli,
#                           ids=high_ids)
# pli_high[[1]]$pli_bin <- 'high'
# pli_high[[2]]$pli_bin <- 'high'
# 
# pli_tot <- rbind.data.frame(pli_run[[1]],
#                             pli_high[[1]])
# 
# pli_combn <- rbind.data.frame(pli_run[[2]],
#                               pli_high[[2]])

# Save 4 regression objects
# total hets, prop AF bins, total hets by pLI, prop AF bins by pLI
save(full_jack_estimate,
     bin_jack_estimate,
     # pli_tot, pli_combn, 
     file=out_obj)
