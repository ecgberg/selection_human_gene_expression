#!/usr/bin/Rscript
require(dplyr)
require(tidyr)

args <- commandArgs(trailingOnly=TRUE)
gnom_file <- args[1]
gtex_file <- args[2]
merged_file <- args[3]

gnomads <- read.table(gnom_file, header=FALSE, sep=' ')
names(gnomads) <- c('CHROM', 'POS', 'REF', 'ALT', 'EUR_AF')
  
gtex <- read.table(gtex_file, sep=' ', skip=1, header=FALSE)
header <- scan(gtex_file,
               nlines=1, what=character())
names(gtex) <- header
  
merged <- left_join(gtex, 
                    gnomads, by=c('#CHROM'='CHROM', 
                                  'POS'='POS',
                                  'REF'='REF',
                                  'ALT'='ALT'))
  
write.table(merged, 
            file=merged_file,
            row.names=FALSE,
            sep=' ',
            quote=FALSE)
