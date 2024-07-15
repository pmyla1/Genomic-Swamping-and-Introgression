#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40g
#SBATCH --time=12:00:00
#SBATCH --job-name=CreateConsensusSequences
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


cd ~/ 

mkdir -p 310524_New_Reference/
###make environmental variables for the reference genome, the input VCF and the F1 output VCF
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/reheadered.F4_133.ann.vcf.gz
DIP=~/310524_New_Reference/310524_UKDiploidsOnly.vcf.gz
TET=~/310524_New_Reference/310524_UKTetraploidsOnly.vcf.gz
HEX=~/310524_New_Reference/310524_UKHexaploidsOnly.vcf.gz

################
##load samtools to index the reference genome 
module load samtools-uoneasy/1.18-GCC-12.3.0

samtools faidx $REF

module unload samtools-uoneasy/1.18-GCC-12.3.0
##############

################
##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##UK DIPLOIDS ONLY
##GATK SelectVariants for UK Diploids only - exclude insertion/deletion mutations, exclude variants with mixed SNPs/indels at the same site
##restrict alleles to biallelic SNPs.
gatk SelectVariants \
 -R $REF \
 -V $VCF \
 -O $DIP \
 -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 \
 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 \
 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 \
 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 \
 --select-type-to-exclude INDEL \
 --select-type-to-exclude MIXED \
 --restrict-alleles-to BIALLELIC \
 --select "AF > 0.25 && AF < 1.0" 
##############   
##UK TETRAPLOIDS ONLY
##GATK SelectVariants for UK TETRAPLOIDS only - exclude insertion/deletion mutations, exclude variants with mixed SNPs/indels at the same site
##restrict alleles to biallelic SNPs.
gatk SelectVariants \
 -R $REF \
 -V $VCF \
 -O $TET \
 -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 \
 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 \
 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 \
 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn GEO_2 -sn GEO_6 \
 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 \
 -sn LOS_1 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 \
 -sn SCU_1 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 \
 -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 \
 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 \
 --select-type-to-exclude INDEL \
 --select-type-to-exclude MIXED \
 --restrict-alleles-to BIALLELIC \
 --allow-nonoverlapping-command-line-samples \
 --select "AF > 0.25 && AF < 1.0" 
#############   
##############   
##UK HEXAPLOIDS ONLY
##GATK SelectVariants for UK HEXAPLOIDS only - exclude insertion/deletion mutations, exclude variants with mixed SNPs/indels at the same site
##restrict alleles to biallelic SNPs.
gatk SelectVariants \
 -R $REF \
 -V $VCF \
 -O $HEX \
 -sn CUM_1 -sn FOR_1 -sn FRE_013 -sn DAR_1 -sn DAR_3 -sn JON_001 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 -sn RYE_1 -sn SCO_1 \
 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 --select-type-to-exclude INDEL \
 --select-type-to-exclude MIXED \
 --restrict-alleles-to BIALLELIC \
 --select "AF > 0.25 && AF < 1.0" 
############## 
###############
##create a sequence dictionary for the reference genome
gatk CreateSequenceDictionary -R $REF 

##unload the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##########

echo "DONE!!"

##lastline
