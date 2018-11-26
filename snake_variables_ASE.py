import os, sys

##########
# RUN VARIABLES
##########

Rscript = "Rscript --no-restore --no-save "
DATA_PATH = '/srv/gsfs0/projects/pritchard/Ziyue/Exp_Sel/GTEx_v6'
TMP_PATH = '/scratch/users/eglassbe'
SLURM_PATH = TMP_PATH + '/slurm_files'
BASE_PATH = '/srv/gsfs0/projects/pritchard/Emily/ASE'
RESULTS_PATH = BASE_PATH + '/results'
SCRIPT_PATH = BASE_PATH + '/smscripts'
PLI_FN = '/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'
DATASET = 'GTExV6'
DATE = '112018'

RACE_CODE = '3' # This is code for self-reported 'white' ancestry from phs000424.v4.pht002742.v4.GTEx_Subject_Phenotypes.data_dict.xml
# grep ETHN on above file (in Ziyue's Sample_Info drive on scg4) for other race/ethnicity codes

TOTAL_READS=5
EXP_TISSUE=2
MIN_DISTANCE=2

HLA_START=25892529
HLA_END=33436144
HLA_BUFFER=2000000
HWE_THRESH=7.88 #corresponds to p-val of 0.005 (0.5% expected filtered sites, on c21 see nearly 3% of sites filtered)

def get_EUR_samples(vcf,
                    ASEinfo,
                    subject_info):
    import gzip
    import subprocess
    import re

    # Pull out samples for which we have WGS data
    print('checking WGS data')
    with gzip.open(vcf, 'rb') as f_vcf:
        # Skip first lines
        for l, line in enumerate(f_vcf):
            if l>140:
                break

    line = line.decode(encoding='UTF-8')
    subject_names = [t[0:9] for t in line.split('\t')[9:]]

    # Pull out samples for which we have RNA data
    print('checking RNASeq files')
    RNASeq=subprocess.Popen(['tar', '-tzf', ASEinfo], stdout=subprocess.PIPE, shell=False)
    RNASeq_files=RNASeq.stdout.read().decode(encoding='UTF-8')
    RNASeq_files=RNASeq_files.strip('\n').split('\n')
    RNASeq_names=[re.search('GTEX-.{4}', sample).group() for sample in RNASeq_files]

    # Pull out Europeans for further analysis
    with open(subject_info, 'r') as f_subj:
        file_subj=f_subj.read()

    print('pulling out Europeans with both WGS+RNASeq')
    array_subj=[ind.split('\t') for ind in file_subj.strip('\n').split('\n')]
    race_info=array_subj[0].index('RACE')
    white_sample_names=[ ind[1] for ind in array_subj if ind[race_info]==RACE_CODE and ind[1] in subject_names and ind[1] in RNASeq_names ]

    return(white_sample_names)


# ASE_SAMPLES = get_EUR_samples(vcf=DATA_PATH + '/WGS/GTEx_Analysis_2015-01-12_WholeGenomeSeq_148Indiv_GATK_HaplotypeCaller.vcf.gz',
#                               ASEinfo=DATA_PATH + '/ASE_data/GTEX_V6_EXOME_ASE.tar.gz',
#                               subject_info=DATA_PATH + '/Subject_Info/GTEx_SUBJECT_ETHN_INFO')
#
# with open(BASE_PATH+'/asesamples.txt', 'w') as f_sampout:
#     ase_samp_string = '\n'.join([str(i) for i in ASE_SAMPLES ])
#     f_sampout.write(ase_samp_string)

with open(BASE_PATH+'/asesamples.txt', 'r') as f_sampin:
          ASE_SAMPLES = [ ase_line.strip('\n') for ase_line in f_sampin.readlines() ]

TISSUES = [ "ADPSBQ",
            "ARTTBL",
            "LUNG",
            "HRTLV",
            "MSCLSK",
            "NERVET",
            "SKINS",
            "THYROID",
            "WHLBLD" ]

