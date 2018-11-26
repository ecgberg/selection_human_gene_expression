# snakemake file for running ASE analysis

# variables file
include: "snake_variables_ASE.py"

wildcard_constraints:
    ase_sample = 'GTEX-.{4}',
    tissue = '[A-Z]{4,7}',
    chrom = '[0-9]{1,2}'

# 
rule all:
    input:
        expand(BASE_PATH + '/' + DATASET + '_regs_{DIST}_' + DATE + '.RDa' , DIST=[10000, 50000]), # For combined-tissue ASE
        expand(BASE_PATH + '/' + DATASET + '_perms_{DIST}_' + DATE + '.RDa', DIST=[10000, 50000]),
        expand(BASE_PATH + '/gnomAD_regs_{DIST}_' + DATE + '.RDa', DIST=[10000, 50000]),
        expand(BASE_PATH + '/gnomAD_perms_{DIST}_' + DATE + '.RDa', DIST=[10000, 50000])

# -------- PROCESS ASE DATA -------- #

# extract relevant individuals from the GTEX ASE tar
rule extract_individual:
    input:
        DATA_PATH + '/ASE_data/GTEX_V6_EXOME_ASE.tar.gz'
    output:
        temp(TMP_PATH + '/{ase_sample}.exome.ase.table.tsv') # make this temp
    params:
        job_name='extract_{ase_sample}',
        run_time='00:30:00',
        cores='1',
        memory='8G',
        error_out_file=SLURM_PATH + '/extract_{ase_sample}',
        queue='standard'
    shell:
        'cd {TMP_PATH}; '
        'tar -zxvf {input} {wildcards.ase_sample}.exome.ase.table.tsv'

# filter potential genotyping errors from ASE data
rule filter_hets:
    input:
        TMP_PATH + '/{ase_sample}.exome.ase.table.tsv'
    output:
        temp(TMP_PATH + '/filter.{ase_sample}.exome.ase.table.tsv')
    params:
        job_name='filter_{ase_sample}',
        run_time='00:10:00',
        cores='1',
        memory='8G',
        # if this folder doesn't exist v bad things happen
        error_out_file=SLURM_PATH + '/filter_{ase_sample}',
        queue='standard'
    shell:
        'module load r; '
        '{Rscript} {SCRIPT_PATH}/filter_geno_error.R {input} {EXP_TISSUE} {TOTAL_READS} {TISSUES} {output}'

# pull apart individual ASE data by tissue
rule sep_tissue:
    input:
        TMP_PATH + '/filter.{ase_sample}.exome.ase.table.tsv'
    output:
        temp(TMP_PATH + '/{tissue}.{ase_sample}.exome.ase.table.tsv')
    params:
        job_name='{tissue}_{ase_sample}',
        run_time='00:10:00',
        cores='1',
        memory='8G',
        error_out_file=SLURM_PATH + '/{tissue}_{ase_sample}',
        queue='standard'
    shell:
        "echo {wildcards.tissue}; "
        "cat {input} | awk '!/^[MX]/" # filter out mito and sex chromosomes
        "{{if (NR==1 || ($8==\"{wildcards.tissue}\"" # then split by current tissue
        "&& $22==0 && $23==0 && $24==0)) print $0}}' > {output}"  # also check for bad genotype quality warnings

# combine all ASE data (across individuals) for a single tissue
rule agg_tissue:
    input:
        expand(TMP_PATH + '/{{tissue}}.{ase_sample}.exome.ase.table.tsv', ase_sample=ASE_SAMPLES)
    output:
        TMP_PATH + '/{tissue}.table.tsv'
    params:
        job_name='{tissue}_table',
        run_time='00:10:00',
        cores='1',
        memory='8G',
        error_out_file=SLURM_PATH + '/{tissue}_table',
        queue='standard'
    shell:
        "cat {input} > {TMP_PATH}/tmp{wildcards.tissue}; "
        "awk '{{if (NR==1 || $1!=\"CHR\") print $0}}' {TMP_PATH}/tmp{wildcards.tissue} > {output}; " #  Remove header lines from each chromosome file
        "rm {TMP_PATH}/tmp{wildcards.tissue} "


# summarize ASE data across tissues to get combined_tissue ASE
## for each gene, selects the best-sampled ASE site (separate job does this for the tissue-specific cases)
rule cross_tissue_ASE:
    input:
        expand(BASE_PATH + '/tissue_level/filterHWE.' + DATASET + '_{tissue}_' + DATE + '.txt', tissue=TISSUES)
    output:
        BASE_PATH + '/full_' + DATASET + '_cross_tissueASE_' + DATE + '.txt'
    params:
        job_name='cross_tissue_ASE',
        run_time='01:00:00',
        cores='1',
        memory='8G',
        error_out_file=SLURM_PATH + '/cross_tissue',
        queue='standard'
    shell:
        'module load r; '
        '{Rscript} {SCRIPT_PATH}/add_tissues.R {input} {TISSUES} {output}'


# separate GTEx WGS vcf (called w/ GATK) by chromosome
rule get_vcf_bychrom:
    input:
        DATA_PATH + '/WGS/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz'
    output:
        TMP_PATH + '/' + DATASET + '_VCF_{chrom}.txt'
    params:
        job_name='c{chrom}_VCF',
        run_time='02:00:00',
        cores='1',
        memory='16G',
        # if this folder doesn't exist v bad things happen
        error_out_file=SLURM_PATH + '/c{chrom}_VCF',
        queue='standard'
    shell:
        "zcat {input} | awk '{{if (NR==142 || $1==\"{wildcards.chrom}\") print $0}}' > {output}" # put header line and relevant genotypes


# Get metadata file with ASE genes
rule get_ASEgenes:
    input:
        ase_sites=BASE_PATH + '/full_' + DATASET + '_cross_tissueASE_' + DATE + '.txt',
        imprint_genes=BASE_PATH + '/human_imprinted_genes.tsv',
        pli_info=BASE_PATH + PLI_FN 
    output:
        BASE_PATH + '/' + DATASET + '_cross_tissueASE_' + DATE + '.txt'
    params:
        job_name='ase_metadata',
        run_time='02:00:00',
        cores='1',
        memory='4G',
        error_out_file=SLURM_PATH + '/asemeta',
        queue='standard'
    shell:
        'module load r; '
        '{Rscript} {SCRIPT_PATH}/get_ASE_gene.R {input.ase_sites} '
        '{HLA_START} {HLA_END} {HLA_BUFFER} {input.imprint_genes} '
        '{input.pli_info} {output}'

# for one chrom, SNPs that are within range of the TSS of a gene that we have ASE info for
# check for HWE
# and parse file into hard genotype calls
rule get_proximal_SNPs:
    input:
        vcf = TMP_PATH + '/' + DATASET + '_VCF_{chrom}.txt',
        asemeta = BASE_PATH + '/' + DATASET + '_cross_tissueASE_' + DATE + '.txt'
    output:
        vcf = TMP_PATH + '/f.' + DATASET + '_VCF_crossTISSUE_{DIST}_{chrom}_' + DATE +'.txt',
        hwe_stats = TMP_PATH + '/HWEstats.' + DATASET + '_{DIST}_{chrom}_' + DATE + '.txt'
    params:
        job_name='filter.c{chrom}_VCF',
        run_time='08:00:00',
        cores='1',
        memory='16G',
        error_out_file=SLURM_PATH + '/filter.c{chrom}_VCF',
        queue='extended'
    run:
        import numpy as np
        from collections import Counter
        with open(input.asemeta) as meta_file:
            gene_table=meta_file.read().strip('\n')
            gene_array=[ line.split(' ') for line in gene_table.split('\n') ]
            # Pull just the current chromosome
            np_gene_array=np.array(gene_array)
            chrom_array=np_gene_array[ np_gene_array[:,5]==wildcards.chrom ] # remove subject ID/ASE columns, get uniques, use after that
            slimmed_genes=np.unique(chrom_array[:,[0,2,4,6]], axis=0)
        with open(output.vcf, 'w') as out_vcf_file, open(output.hwe_stats, 'w') as hwe_stats:
            with open(input.vcf) as vcf_file:
                # Add appropriate headers from metadata and VCFs
                header_line=vcf_file.readline().strip('\n')
                vcf_head=header_line.split('\t')
                # get indices of European samples to calculate allele frequencies
                ## removing the metadata columns, as they are not included in the genotype calculations
                eur_indices = [i for i, samp_id in enumerate(vcf_head) if any(sample in samp_id for sample in ASE_SAMPLES)]
                eur_sample_ids = [ vcf_head[i] for i in eur_indices ]
                header=' '.join([str(i) for i in np_gene_array[0][[0,2,4,6]].tolist() + vcf_head[0:5] + ['EUR_ALT_COUNTS', 'EUR_MIN_COUNTS', 'EUR_SAMPLE_SIZE'] + eur_sample_ids ])+'\n'
                out_vcf_file.write(header)
                stat_header=' '.join(vcf_head[0:2]) + ' '+ ' '.join([str(i) for i in np_gene_array[0]]) + ' DIST_TO_TSS FREQ CHISQ\n'
                hwe_stats.write(stat_header)
                for line in vcf_file:
                    line=line.strip('\n')
                    SNP=line.split('\t')
                    if (len(SNP)>1): # check for empty line at EOF
                        for ASE_gene in slimmed_genes:
                            # calculate distance from SNP to transcription start site
                            dist = abs( int(SNP[1]) - int(ASE_gene[3]) )
                            # if this SNP is within the correct distance AND
                            # if this SNP is a true SNP (not an indel), parse line
                            if (dist < int(wildcards.DIST) and dist > MIN_DISTANCE
                                and len(SNP[3])==1 and len(SNP[4])==1):
                                EUR_genotypes = [ SNP[i] for i in eur_indices ] # Pull out ONLY genotypes of Europeans
                                genotypes = [ genotype.split(':')[0] for genotype in EUR_genotypes ]
                                geno_calls = [ geno.split('/') for geno in genotypes ]
                                numeric_genotypes = [ int(g[0]) + int(g[1]) if (g[0]!='.' or g[1]!='.') else 'nan' for g in geno_calls ]
                                np_genotypes = np.asarray( numeric_genotypes, dtype=float )
                                eur_sample_size = sum([ (not np.isnan(geno)) for geno in np_genotypes ])*2
                                eur_altcounts = int(np.nansum( np_genotypes ))
                                eur_mincounts = min( eur_altcounts, (eur_sample_size - eur_altcounts) )
                                # Test for HWE
                                obs=Counter(numeric_genotypes)
                                print(obs)
                                obs.pop('nan', None)
                                if (sum(obs.values())>0):
                                    freq=float(sum(filter(lambda i: isinstance(i, int), numeric_genotypes))) / (sum(obs.values())*2)
                                else:
                                    continue
                                    print('done')
                                exp={0:(1-freq)**2 * sum(obs.values()), 1:2*freq*(1-freq) * sum(obs.values()), 2:freq**2 * sum(obs.values())}
                                chisq=sum( [ (obs[i]-exp[i])**2 / exp[i] if exp[i]>0 else 0 for i in [0,1,2] ] )
                                print(chisq)
                                hwe_stats.write(' '.join([str(i) for i in SNP[0:2]]) + ' ' + ' '.join([str(i) for i in ASE_gene]) + ' ' + str(dist) + ' ' + str(freq) + ' ' + str(chisq) + '\n')
                                # If it is in HWE, write all individuals with that variant to the filtered file
                                # If this position is variable within GTEx in Europe and in HWE, write to a file
                                if (eur_mincounts>0 and chisq<HWE_THRESH):
                                    newline=' '.join([str(i) for i in ASE_gene.tolist() + SNP[0:5] + [ eur_altcounts, eur_mincounts, eur_sample_size ] + numeric_genotypes])+'\n'
                                    out_vcf_file.write(newline)


# for one chrom at a time, filter ASE sites by HWE
rule filter_HWE:
    input:
        vcf_files = expand(TMP_PATH + '/' + DATASET + '_VCF_{chrom}.txt', chrom=range(1,23)),
        ASE_file = TMP_PATH + '/{tissue}.table.tsv' 
    output:
        filtered_sites = BASE_PATH + '/tissue_level/filterHWE.' + DATASET + '_{tissue}_' + DATE + '.txt',
        stat_file = TMP_PATH + '/HWE_stats_{tissue}' + DATE + '.txt'
    params:
        job_name='HWE_{tissue}_VCF',
        run_time='6:00:00',
        cores='1',
        memory='16G',
        error_out_file=SLURM_PATH + '/HWE_{tissue}_VCF',
        queue='standard'
    run:
        from collections import Counter
        with open(output.filtered_sites, 'w') as out_sitefile, open(output.stat_file, 'w') as out_statfile:
            with open(input.ASE_file) as ASE_file:
                # Set up file for filtered ASE sites
                header=ASE_file.readline()
                out_sitefile.write(header)
                # And one to write down the chisq stats for each of the sites tested
                stat_header='SITE FREQ OBS_HOM1 OBS_HET OBS_HOM2 CHISQ\n'
                out_statfile.write(stat_header)
                # Build a dictionary of ASE sites to check in the VCF file
                ASE_sites = dict()
                for site in ASE_file:
                    split_site=site.strip('\n').split(' ')
                    site_key=':'.join(split_site[0:3])
                    try:
                        add_line=ASE_sites[site_key]+' '.join(split_site)+'\n'
                        ASE_sites[site_key]=add_line # If that chrom+position are already in the dict, add the line
                    except KeyError:
                        ASE_sites[site_key]=' '.join(split_site)+'\n'
                for vcf_file in input.vcf_files:
                    with open(vcf_file) as VCF_file:
                        # get indices of European samples to calculate allele frequencies
                        ## removing the metadata columns, as they are not included in the genotype calculations
                        header_line=VCF_file.readline()
                        vcf_head=header_line.split('\t')
                        eur_indices = [i for i, samp_id in enumerate(vcf_head) if any(sample in samp_id for sample in ASE_SAMPLES)]
                        # Loop over VCF, only write out ASE sites that are in HWE
                        for line in VCF_file:
                            line=line.strip('\n').split('\t')
                            curr_site=':'.join(line[0:3])
                            try:
                                ASE_sites[curr_site] # Check if this site is worth checking at this point
                                EUR_genotypes = [ line[i] for i in eur_indices ]
                                genotypes = [ genotype.split(':')[0] for genotype in EUR_genotypes ]
                                geno_calls = [ geno.split('/') for geno in genotypes ]
                                numeric_genotypes = [ int(g[0]) + int(g[1]) if (g[0]!='.' or g[1]!='.') else 'nan' for g in geno_calls ]
                                obs=Counter(numeric_genotypes)
                                print(obs)
                                obs.pop('nan', None)
                                if (sum(obs.values())>0):
                                    freq=float(sum(filter(lambda i: isinstance(i, int), numeric_genotypes))) / (sum(obs.values())*2)
                                else:
                                    continue
                                    print('done')
                                exp={0:(1-freq)**2 * sum(obs.values()), 1:2*freq*(1-freq) * sum(obs.values()), 2:freq**2 * sum(obs.values())}
                                chisq=sum( [ (obs[i]-exp[i])**2 / exp[i] if exp[i]>0 else 0 for i in [0,1,2] ] )
                                print(chisq)
                                # If it is in HWE, write all individuals with that variant to the filtered file
                                if(chisq<HWE_THRESH):
                                    out_sitefile.write( ASE_sites[curr_site] )
                                # Regardless, write that stat to a stat file
                                stat_line=curr_site + ' ' + ' '.join([ str(val) for val in [freq, obs[0], obs[1], obs[2], chisq] ]) +'\n'
                                out_statfile.write(stat_line)
                            except KeyError:
                                continue
         

# Figure out appropriate allele frequency bin breaks for regression
## i.e. count number of singletons and allele frequency bins
rule count_bins:
    input:
        vcfs = expand(TMP_PATH + '/GTExV6_VCF_{{DIST}}_{chrom}_gnomAD.txt', chrom=range(1,23))
    output:
        gtex = BASE_PATH + '/gtex_bins_{DIST}_' + DATE,
        gnomad = BASE_PATH + '/gnomad_bins_{DIST}_' + DATE
    params:
        job_name='count_bin_breaks',
        run_time='03:00:00',
        cores='1',
        memory='36G',
        # if this folder doesn't exist v bad things happen
        error_out_file=SLURM_PATH + '/count_bin_{wildcards.dist}.breaks',
        queue='standard'
    run:
        import csv
        import numpy as np
        import math
        from collections import Counter
        from heapq import nsmallest

        # Loop over each file, pull out GTEx minor allele counts and gnomAD frequencies
        minor_allele_counts = []
        gtex_maf = []
        gnom_maf = []
        for vcf in input.vcfs:

            with open(vcf) as c_vcf:
                read_vcf = csv.reader(c_vcf, delimiter=' ')
                t = list(zip(*read_vcf))

            # Store this chromosome's minor allele information
            c_counts = [ int(counts) for counts in t[10][1:] ]
            c_samp = [ float(samp) for samp in t[11][1:] ] # Need to track sample size bc some individuals aren't confidently genotyped at every position
            # Calculate GTEx minor allele frequencies (for binning above doubletons)
            c_maf = [ count/samp for count, samp in zip(c_counts, c_samp) if count > 2 ]
            c_gnom = [min(float(af), 1-float(af)) for af in t[134][1:] if (af!='NA' and af!='.')]

            minor_allele_counts = minor_allele_counts + c_counts
            gtex_maf = gtex_maf + c_maf
            gnom_maf = gnom_maf + c_gnom

            print('done with '+vcf)

        # COUNT GTEx BINS
        MAc = Counter(minor_allele_counts)
        doubletons = MAc[2]

        # Initialize a vector to store bin breaks based on MAF counts
        # There are 122 European samples, so that's the max number of minor alleles
        # Initialize with a singleton bin
        i=1
        gtex_breaks=[0]
        while gtex_breaks[-1]<0.5:
            gtex_breaks = gtex_breaks + [ nsmallest(doubletons*i, gtex_maf)[-1] ]
            i+=1
        
        with open(output.gtex, 'w') as outgtex:
            outgtex.write(' '.join( [str(b) for b in gtex_breaks ]))

        gnom_breaks=np.percentile(gnom_maf, np.arange(0,110,11.111111))
        
        with open(output.gnomad, 'w') as outgnomad:
            outgnomad.write(' '.join( [str(b) for b in gnom_breaks ]))



# Count the number of cis-heterozygous sites in each allele frequency bin around each ASE site
## Returns a file ready for regression (column of adjusted J2)
## Counts spectra binned by gnomAD and by GTEx
rule count_spectra:
    input:
        vcf = TMP_PATH + '/' + DATASET + '_VCF_{DIST}_{chrom}_gnomAD.txt', 
        ase_sites = BASE_PATH + '/' + DATASET + '_cross_tissueASE_' + DATE + '.txt',
        gtex_bins = BASE_PATH + '/gtex_bins_{DIST}_' + DATE,
        gnomAD_bins = BASE_PATH + '/gnomad_bins_{DIST}_'+DATE # text file with bin breaks
    output:
        unique_variants=TMP_PATH + '/' + DATASET + 'uniques_crossTISSUE_{chrom}_{DIST}' + DATE + '.RDa',
        variants_pergene=TMP_PATH + '/' + DATASET + 'pergenes_crossTISSUE_{chrom}_{DIST}_' + DATE + '.txt',
        gtex_bins=temp(TMP_PATH + '/' + DATASET + '_AFbins_{chrom}_{DIST}_' + DATE + '.txt'), # for summary across tissues
        gnomAD_bins=temp(TMP_PATH + '/gnomAD_AFbins_{chrom}_{DIST}_' + DATE + '.txt') # summary across tissues, binned by gnomAD frequencies
    params:
        job_name='count_c{chrom}_{DIST}spectra',
        run_time='06:00:00',
        cores='1',
        memory='24G',
        # if this folder doesn't exist v bad things happen
        error_out_file=SLURM_PATH + '/count_{chrom}_spectra',
        queue='standard'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/count_AFspectra.R {input.vcf} {input.ase_sites} {wildcards.chrom} {input.gtex_bins} {input.gnomAD_bins} {output.unique_variants} {output.variants_pergene} {output.gtex_bins} {output.gnomAD_bins}"

# # put AF spectrum counts for all ASE-containing individual-gene pairs
rule merge_spectra:
    input:
        gtex_chr_spectra=expand(TMP_PATH + '/' + DATASET + '_AFbins_{chrom}_{{DIST}}_' + DATE + '.txt', chrom=range(1,23)),
        gtex_pergene_variants=expand(TMP_PATH + '/' + DATASET + 'pergenes_crossTISSUE_{chrom}_{{DIST}}_' + DATE + '.txt', chrom=range(1,23)),
        gnomAD_chr_spectra=expand(TMP_PATH + '/gnomAD_AFbins_{chrom}_{{DIST}}_' + DATE + '.txt', chrom=range(1,23))
    output:
        gtex_spectra=BASE_PATH + '/' + DATASET + '_AFbins_{DIST}_' + DATE + '.txt', # summarized across tissues
        gnomAD_spectra=BASE_PATH + '/gnomAD_AFbins_{DIST}_' + DATE + '.txt'
    params:
        job_name='merge_spectra',
        run_time='01:00:00',
        cores='1',
        memory='16G',
        # if this folder doesn't exist v bad things happen
        error_out_file=SLURM_PATH+ '/merge_spectra',
        queue='standard'
    shell:
        "cat {input.gtex_chr_spectra} > {TMP_PATH}/tmpGTEx{wildcards.DIST}spectra; " 
        "cat {input.gnomAD_chr_spectra} > {TMP_PATH}/tmpgnomAD{wildcards.DIST}spectra; " 
        "awk '{{if (NR==1 || $1 !~ \"SUBJECT_ID\") print $0}}' {TMP_PATH}/tmpGTEx{wildcards.DIST}spectra > {output.gtex_spectra}; " 
        "awk '{{if (NR==1 || $1 !~ \"SUBJECT_ID\") print $0}}' {TMP_PATH}/tmpgnomAD{wildcards.DIST}spectra > {output.gnomAD_spectra}; " 
        "rm {TMP_PATH}/tmp*{wildcards.DIST}spectra"

# perform AF bin regressions - overall AND by pLI bin - jackknife confidence intervals shown in Figure 6
# For GTEx allele frequency bins
## including jackknife to estimate confidence intervals
rule run_gtex_ASEregs:
    input:
        gtex_AF_bincounts=BASE_PATH + '/' + DATASET + '_AFbins_{DIST}_' + DATE + '.txt',
        # AF_bincounts=BASE_PATH + '/ASE/' + DATASET + '_AFbins_{tissue}_{DIST}_' + DATE + '.txt',
    output:
        gtex_regs=BASE_PATH + '/' + DATASET + '_regs_{DIST}_' + DATE + '.RDa' # will be an object in 4 pcs (total, total by pLI, by AF bin, AF bin by pLI)
        # BASE_PATH + '/ASE/' + DATASET + '_{tissue}regs_{DIST}_' + DATE + '.RDa' # will be an object in 4 pcs (total, total by pLI, by AF bin, AF bin by pLI)
    params:
        job_name='regress_tissue',
        run_time='8:00:00',
        cores='1',
        memory='12G',
        error_out_file=SLURM_PATH+ '/regress_tissue',
        queue='extended'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/jack_reg.R {input.gtex_AF_bincounts} {output.gtex_regs}" # run for gtex


# perform AF bin regressions - overall AND by pLI bin - jackknife confidence intervals shown in Figure 6
# For GTEx allele frequency bins
## including PERMUTATION TEST
rule run_gtex_ASEperms:
    input:
        gtex_AF_bincounts=BASE_PATH + '/' + DATASET + '_AFbins_{DIST}_' + DATE + '.txt',
        # AF_bincounts=BASE_PATH + '/ASE/' + DATASET + '_AFbins_{tissue}_{DIST}_' + DATE + '.txt',
    output:
        gtex_perms=BASE_PATH + '/' + DATASET + '_perms_{DIST}_' + DATE + '.RDa' # will be an object in 4 pcs (total, total by pLI, by AF bin, AF bin by pLI)
        # BASE_PATH + '/ASE/' + DATASET + '_{tissue}regs_{DIST}_' + DATE + '.RDa' # will be an object in 4 pcs (total, total by pLI, by AF bin, AF bin by pLI)
    params:
        job_name='permute_tissue',
        run_time='8:00:00',
        cores='1',
        memory='12G',
        error_out_file=SLURM_PATH + '/permute_tissue',
        queue='extended'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/permute_regs.R {input.gtex_AF_bincounts} {output.gtex_perms}" # run for gtex

# perform AF bin regressions - overall AND by pLI bin - jackknife confidence intervals shown in Figure 6
# For gnomAD allele frequency bins
## including jackknife to estimate confidence intervals
rule run_gnom_ASEregs:
    input:
        gnomAD_spectra=BASE_PATH + '/gnomAD_AFbins_{DIST}_' + DATE + '.txt',
        # AF_bincounts=BASE_PATH + '/ASE/' + DATASET + '_AFbins_{tissue}_{DIST}_' + DATE + '.txt',
    output:
        gnomAD_regs=BASE_PATH + '/gnomAD_regs_{DIST}_' + DATE + '.RDa'
    params:
        job_name='regress_gnom_tissue',
        run_time='8:00:00',
        cores='1',
        memory='12G',
        error_out_file=SLURM_PATH + '/regress_gnom_tissue',
        queue='extended'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/jack_reg_gnom.R {input.gnomAD_spectra} {output.gnomAD_regs}" # run regressions for gnomAD data


# perform AF bin regressions - overall AND by pLI bin - jackknife confidence intervals shown in Figure 6
# For gnomAD allele frequency bins
## including permutation test
rule run_gnom_ASEperms:
    input:
        gnomAD_spectra=BASE_PATH + '/gnomAD_AFbins_{DIST}_' + DATE + '.txt',
        # AF_bincounts=BASE_PATH + '/ASE/' + DATASET + '_AFbins_{tissue}_{DIST}_' + DATE + '.txt',
    output:
        gnomAD_perms=BASE_PATH + '/gnomAD_perms_{DIST}_' + DATE + '.RDa'
    params:
        job_name='permute_gnom_tissue',
        run_time='8:00:00',
        cores='1',
        memory='12G',
        error_out_file=SLURM_PATH + '/permute_gnom_tissue',
        queue='extended'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/permute_reg_gnom.R {input.gnomAD_spectra} {output.gnomAD_perms}" # run regressions for gnomAD data



##### ------ PROCESS gnomAD allele frequencies -------- #####

# For each chrom, list variants that appear in the GTEx sample
## To query these allele frequencies in gnomAD
rule variants_from_GTEx:
    input:
        TMP_PATH + '/f.' + DATASET + '_VCF_crossTISSUE_{DIST}_{chrom}_' + DATE +'.txt'
    output:
        temp(TMP_PATH + '/' + DATASET + '_variants_{DIST}_{chrom}.tsv')
    params:
        job_name='GTEx_variants',
        run_time='08:00:00',
        cores='1',
        memory='12G',
        error_out_file=SLURM_PATH+'/list_variants',
        queue='standard'
    shell:
        "awk '{{ if (NR!=1) print $5\"\t\"$6}}' {input} > {output}"

# For each chrom, order the variants from the GTEx sample
## Generate region file to be used with BCFtools
rule order_GTEx_variants:
    input:
        TMP_PATH + '/' + DATASET + '_variants_{DIST}_{chrom}.tsv'
    output:
        temp(TMP_PATH + '/' + DATASET + '_ordered_variants_{DIST}_{chrom}.tsv')
    params:
        job_name='order_GTEx_variants',
        run_time='04:00:00',
        cores='1',
        memory='6G',
        error_out_file=SLURM_PATH + '/order_variants{wildcards.chrom}',
        queue='standard'
    shell:
        "module load r; "
        "{Rscript} {SCRIPT_PATH}/order_variants.R {input} {output}"

# For each chrom, pull out the gnomAD European allele frequencies
## Each ref/alt allele on one line
rule get_gnomAD_freqs:
    input:
        gtex_variants=TMP_PATH + '/' + DATASET + '_ordered_variants_{DIST}_{chrom}.tsv',
        gnomAD_VCF=BASE_PATH + '/gnomADs/gnomad.genomes.r2.0.2.sites.chr{chrom}.vcf.bgz'
    output:
        TMP_PATH + '/gnomAD_AFs_{DIST}_{chrom}.tsv'
    params:
        job_name='get_gnomAD_AFs',
        run_time='01:00:00',
        cores='1',
        memory='6G',
        error_out_file=SLURM_PATH + '/gnomAD_afs',
        queue='standard'
    shell:
        '{BASE_PATH}/bin/bcftools view -T {input.gtex_variants} -v snps {input.gnomAD_VCF} | {BASE_PATH}/bin/bcftools norm -m -any | {BASE_PATH}/bin/bcftools query -f "%CHROM %POS %REF %ALT %INFO/AF_NFE\n" -o {output}'

# Merge gnomAD frequencies to GTEx (include REF/ALT check)
rule merge_gnomAD:
    input:
        gnomAD_frequencies=TMP_PATH + '/gnomAD_AFs_{DIST}_{chrom}.tsv',
        GTEx_data=TMP_PATH + '/f.' + DATASET + '_VCF_crossTISSUE_{DIST}_{chrom}_' + DATE +'.txt'
    output:
        TMP_PATH + '/' + DATASET + '_VCF_{DIST}_{chrom}_gnomAD.txt' # This is the file we will use to count allele frequency
    params:
        job_name='add_gnomAD_toGTEx',
        run_time='02:00:00',
        cores='1',
        memory='16G',
        error_out_file=SLURM_PATH + '/gnomAD_GTEx',
        queue='standard'
    shell:
        'module load r; '
        '{Rscript} {SCRIPT_PATH}/merge_gnomAD_GTEx.R '
        '{input.gnomAD_frequencies} {input.GTEx_data} {output}'

