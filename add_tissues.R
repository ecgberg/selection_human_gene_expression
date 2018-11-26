#!/usr/bin/Rscript
require(dplyr)
require(foreach)

args <- commandArgs(trailingOnly=TRUE)
list_files <- args[1:(length(args)-10)]
tissues <- args[10:(length(args)-1)]
outfile <- args[length(args)]

individuals <- foreach(fn=list_files, .combine=rbind.data.frame) %do% {
    ind <- read.table(fn, header=TRUE, stringsAsFactors=FALSE)
}

print('data loaded')
summ_tissue <- individuals %>% filter(TISSUE_ID %in% tissues) %>%
                               group_by(CHR, 
                                        POS, 
                                        VARIANT_ID,
                                        VARIANT_ANNOTATION,
                                        REF_ALLELE,
                                        ALT_ALLELE,
                                        SUBJECT_ID,
                                        GENE_ID) %>%
                                summarize(REF_COUNT=sum(REF_COUNT),
                                          ALT_COUNT=sum(ALT_COUNT),
                                          TOTAL_COUNT=sum(TOTAL_COUNT)) %>% 
                                mutate(TISSUE_ID='CROSS-TISSUE')

write.table(summ_tissue, file=outfile, row.names=FALSE, quote=FALSE)
