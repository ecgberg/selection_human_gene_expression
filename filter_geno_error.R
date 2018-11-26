#!/usr/bin/Rscript
## pre-process GTEx ASE table to filter out potential genotyping errors
require(dplyr)

args=commandArgs(trailingOnly=TRUE)
infile <- args[1]
exp_tissue <- as.numeric(args[2])
total_reads <- as.numeric(args[3])
tissues <- args[4:(length(args)-1)]
outfile <- args[length(args)]


ASE_file <- read.table(infile, header=TRUE)

# Filter to ASE files we're interested in
# ASE_tissues <- ASE_file %>% filter(TISSUE_ID %in% tissues)

## For every variant, 
### calculate the number of tissues each allele appears in and 
### the total number of reads (across tissues) supporting that read
filtered_data <- ASE_file %>%
                    group_by(VARIANT_ID) %>%
                    mutate(ref_total=sum(REF_COUNT),
                           alt_total=sum(ALT_COUNT),
                           n_ref_tissue=sum(REF_COUNT>0),
                           n_alt_tissue=sum(ALT_COUNT>0)) %>%
                    filter(ref_total>=total_reads,
                           alt_total>=total_reads,
                           n_ref_tissue>=exp_tissue,
                           n_alt_tissue>=exp_tissue) %>%
                     select(-ref_total, -alt_total, -n_ref_tissue, -n_alt_tissue) %>%
                     ungroup()
                      
# Write to appropriate file
write.table(filtered_data, file=outfile, row.names = FALSE, quote=FALSE)

