#!/usr/bin/rscript
require(dplyr)
require(tidyr)

args <- commandArgs(trailingOnly=TRUE)
vcf_file <- args[1]
tissue_file <- args[2]
chrom <- as.numeric(args[3])
gtex_breaks <- scan(args[4], what=character())
gnomAD_breaks <- scan(args[5], what=character())
gtex_unique_variants <- args[6]
gtex_pergene_hets <- args[7]
gtex_regressionfile <- args[8]
gnomAD_regressionfile <- args[9]

print(args)


## ---------- READ IN AND CLEAN UP BOTH GENOTYPE AND ASE DATA ----------- ##
# -- GENOTYPES FIRST
# Read in header and file, keeping only subject (not sample) ID for columns
header <- scan(vcf_file, what=character(), nlines=1)
header[13:(length(header)-1)] <- substr(header[13:(length(header)-1)], start=1, stop=9)
vcf <- read.table(vcf_file, 
                  sep=' ', header=FALSE, 
                  stringsAsFactors=FALSE, skip=1, na.strings='nan')
names(vcf) <- header
print('VCF file is read in')

# -- AF BIN COUNTING -- gnomAD data
vcf$gnom_EUR_MAF <- sapply(as.numeric(vcf$EUR_AF), function(af) min(af, 1-af))
vcf$gnom_bins <- .bincode(vcf$gnom_EUR_MAF, gnomAD_breaks, include.lowest = TRUE, right=TRUE)

# -- AF BIN COUNTING -- GTEx data
vcf$gtex_AF <- vcf$EUR_MIN_COUNTS/vcf$EUR_SAMPLE_SIZE
vcf$gtex_EUR_MAF=sapply(vcf$gtex_AF, function(af) min(af, 1-af))
vcf$gtex_bins=.bincode(vcf$gtex_EUR_MAF, gtex_breaks, include.lowest=TRUE, right=TRUE) + 2
vcf <- vcf %>% mutate(gtex_bins=ifelse(EUR_MIN_COUNTS<=2, EUR_MIN_COUNTS, gtex_bins))

# -- Count spectra
sample_names <- names(vcf)[grepl('GTEX-', names(vcf))]
gtex_gene_counts <- lapply(sample_names, function(subject) { 
    vcf %>% # Loop over all individuals
        filter(get(subject)==1) %>% # For that individual, pull out all heterozygous sites
        group_by(ensembl_gene_id, 
                 transcription_start_site, 
                 `#CHROM`, gtex_bins) %>% # Count allele frequency bins of SNPs around that TSS
                 tally() %>% 
                 mutate(SUB=substr(subject, 1,9)) 
    } ) %>% # Add a subject identifier for later
    bind_rows %>% # Put rows together for each individual
    spread(gtex_bins, n, fill=0) # Make this in wide format, casting NAs (empty AF bins) to 0s for regression

gnom_gene_counts <- lapply(sample_names, function(subject) { 
    vcf %>% # Loop over all individuals
        filter(get(subject)==1) %>% # For that individual, pull out all heterozygous sites
        group_by(ensembl_gene_id, 
                 transcription_start_site, 
                 `#CHROM`, gnom_bins) %>% # Count allele frequency bins of SNPs around that TSS
                 tally() %>% 
                 mutate(SUB=substr(subject, 1,9)) 
    } ) %>% # Add a subject identifier for later
    bind_rows %>% # Put rows together for each individual
    spread(gnom_bins, n, fill=0) # Make this in wide format, casting NAs (empty AF bins) to 0s for regression


# -- Count number of unique variants around start sites of genes (and what their allele frequencies are)
#    Add the average per individual?
unique_vars <- vcf %>% group_by(`#CHROM`, gtex_bins) %>% 
                       summarize(n_vars=length(unique(POS)))
                     
save(unique_vars,
     file=gtex_unique_variants)

# -- Count heterozygosity etc per gene
gtex_variant_counts <- vcf %>%
                       group_by(`#CHROM`,
                                ensembl_gene_id,
                                gtex_bins) %>%
                       summarize(n_tot=length(POS), # Don't require this, really
                                 AF_factor=sum(2*gtex_EUR_MAF*(1-gtex_EUR_MAF)))
write.table(gtex_variant_counts,
            row.names=FALSE,
            quote=FALSE,
            file=gtex_pergene_hets)


# -- NOW ASE
ASE_dat <- read.table(tissue_file, 
                      sep=' ', header=TRUE, stringsAsFactors=FALSE)
# calculate overall, tissue-level reference bias (should I do this for only relevant chromosome? or only the ASE sites with the best reads? I think broader is better)
ref_bias <- sum(ASE_dat$ALT_COUNT) / sum(ASE_dat$TOTAL_COUNT)
print(ref_bias)

# filter to relevant chromosome
c.ASE_dat <- ASE_dat %>% filter(CHR==chrom)

## Add measurement of ASE
j.ASE_dat <- c.ASE_dat %>% 
                 mutate(J2= ( ALT_COUNT - ref_bias*TOTAL_COUNT )**2 / ( ref_bias*(1-ref_bias)*TOTAL_COUNT ),
                        adj_J2= (J2-1) / (4*(TOTAL_COUNT-1)) ) %>% 
                        dplyr::select(SUBJECT_ID,
                                      CHR, 
                                      ensembl_gene_id,
                                      transcription_start_site, 
                                      strand, 
                                      pli_bin, 
                                      POS, 
                                      VARIANT_ID, 
                                      REF_COUNT, 
                                      ALT_COUNT, 
                                      J2, 
                                      adj_J2) # Removed #TISSUE_ID column for joint analysis

## Merge ASE measurements and variant counts - this is where you'll lose a lot of SNP counts
GTEx_spectra <- left_join(j.ASE_dat, gtex_gene_counts, 
                          by=c('CHR'='#CHROM', 'ensembl_gene_id', 'SUBJECT_ID'='SUB', 'transcription_start_site')) # tss

gnomAD_spectra <- left_join(j.ASE_dat, gnom_gene_counts, 
                            by=c('CHR'='#CHROM', 'ensembl_gene_id', 'SUBJECT_ID'='SUB', 'transcription_start_site')) # tss

## write it out
write.table(GTEx_spectra, file=gtex_regressionfile, row.names=FALSE, sep=' ', quote=FALSE)
write.table(gnomAD_spectra, file=gnomAD_regressionfile, row.names=FALSE, sep=' ', quote=FALSE)
