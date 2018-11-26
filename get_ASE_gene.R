#!/usr/bin/Rscript
require(dplyr)
require(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
ase_gene_fn <- args[1]
HLA_start <- as.numeric(args[2])
HLA_end <- as.numeric(args[3])
HLA_buffer <- as.numeric(args[4])
imprint_file <- args[5]
pli_filename <- args[6]
asegene_output <- args[7]

ase_list <- read.table(ase_gene_fn, sep=' ', header=TRUE, stringsAsFactors=FALSE)

# Remove sites that might cause allelic imbalance by damaging the protein
ase_sites <- ase_list %>% filter(VARIANT_ANNOTATION %in% c('3_prime_UTR_variant',
                                                           'synonymous_variant',
                                                           '5_prime_UTR_variant',
                                                           'stop_retained_variant',
                                                           'intron_variant'))

print(paste('Looking for', length(unique(ase_sites$GENE_ID)), 'gene info'))
# Read in raw gene metadata (from biomaRt + ExAC server, Kaitlin Samocha's work)
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="grch37.ensembl.org",
                  dataset = "hsapiens_gene_ensembl")

# Pull down ONLY protein coding genes on chromosomes 1-22
# and also genes for which we have phenotype info (i.e. readcounts)
protein_IDs <- getBM(attributes=c('hgnc_symbol',
                                  'ensembl_gene_id',
                                  'transcription_start_site',
                                  'chromosome_name',
                                  'strand',
                                  'start_position',
                                  'end_position',
                                  'gene_biotype'),
                      filters=c("biotype",
                                "transcript_biotype",
                                "chromosome_name",
                                "ensembl_gene_id"),
                      values=list('protein_coding',
                                  'protein_coding',
                                  c(1:22),
                                  unique(ase_sites$GENE_ID)),
                      mart=ensembl)

protein_IDs[protein_IDs$hgnc_symbol=='','hgnc_symbol'] <- NA
print('got protein_IDs')

m.protein_IDs <- left_join(protein_IDs, ase_sites, by=c('ensembl_gene_id'='GENE_ID'))
print(paste('ASE sites/TSS: ', dim(m.protein_IDs)))

# Remove all ASE site/TSS pairs for which the TSS is downstream of the ASE site
## Note that this may remove some ASE sites from consideration
f.ASEtss_pairs <- m.protein_IDs %>% 
                      filter(ifelse(strand==1, 
                                    transcription_start_site < POS,
                                    transcription_start_site > POS)) %>%
                      ungroup() %>% 
                      group_by(ensembl_gene_id, 
                               hgnc_symbol,
                               chromosome_name, 
                               SUBJECT_ID,
                               TISSUE_ID) %>% # Then, for each individual-gene pair, pick the ase site with the highest coverage
                      filter(TOTAL_COUNT==max(TOTAL_COUNT)) %>% 
                      ungroup() %>% # Then, for each chosen site, pick the most upstream TSS
                      group_by(TISSUE_ID,
                               SUBJECT_ID,
                               ensembl_gene_id,
                               hgnc_symbol,
                               chromosome_name,
                               POS) %>% 
                      filter(ifelse(strand==1,
                                    transcription_start_site==min(transcription_start_site),
                                    transcription_start_site==max(transcription_start_site))) %>% 
                      ungroup() %>% # Then, remove duplicate ASE sites that might result from the same coverage (pick one randomly)
                      group_by(TISSUE_ID,
                               SUBJECT_ID,
                               ensembl_gene_id,
                               hgnc_symbol,
                               chromosome_name) %>%
                      sample_n(size=1)

print(paste('ASE-site/individual/TSS pairs sampled'))

# Add pLI info
pLI <- read.table(pli_filename,
                  header=TRUE)

### Set up appropriate pLI bins
pLI$pli_bin <- rep('med', nrow(pLI))
pLI[pLI$pLI>0.9, 'pli_bin'] <- 'high'
pLI[pLI$pLI<=0.1, 'pli_bin'] <- 'low'
pLI$pli_bin <- factor(pLI$pli_bin)
pLI$pli_bin <- factor(pLI$pli_bin,
                      levels=levels(pLI$pli_bin)[c(2,3,1)])

pLI <- subset(pLI, select=c('gene', 'pli_bin'))
names(pLI) <- c('hgnc_symbol', 'pli_bin')

# put pLI and TSS info together
gene_metadata <- left_join(f.ASEtss_pairs %>% ungroup() %>% 
                                              dplyr::select(TISSUE_ID,
                                                            SUBJECT_ID,
                                                            ensembl_gene_id,
                                                            hgnc_symbol,
                                                            strand,
                                                            CHR,
                                                            transcription_start_site,
                                                            POS,
                                                            VARIANT_ID,
                                                            REF_ALLELE,
                                                            ALT_ALLELE,
                                                            REF_COUNT,
                                                            ALT_COUNT,
                                                            TOTAL_COUNT), 
                           pLI, 
                           by='hgnc_symbol', 
                           all.x=TRUE)

# REMOVE IMPRINTED AND HLA GENES
imprint_genes <- read.table(imprint_file, header=TRUE, stringsAsFactors = FALSE)
filt_spec <- gene_metadata %>% 
                 filter(!( hgnc_symbol %in% imprint_genes$Gene )) %>%
                 filter(!( CHR==6 & (transcription_start_site > (HLA_start - HLA_buffer) & 
                                       transcription_start_site < (HLA_end + HLA_buffer)) ))

# save files
write.table(filt_spec, file=asegene_output,
            sep=' ', quote=FALSE, row.names=FALSE)
