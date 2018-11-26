#!/usr/bin/Rscript
require(dplyr)
require(tidyr)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

input_variants <- read.table(input_file, header=FALSE, stringsAsFactors=FALSE)
names(input_variants) <- c('CHROM', 'POS')

ord_file <- input_variants[order(input_variants$CHROM,
                                 input_variants$POS),]
ord_file <- ord_file %>% distinct()

write.table(ord_file, file=output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
