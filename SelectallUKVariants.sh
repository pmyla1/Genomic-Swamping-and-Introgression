#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=gatk_select_allUKtets_allUKhex
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

####################
##This script uses GATK SelectVariants (Version 4.4.0) to select all UK hexaploids, including Cochlearia danica and Cochlearia anglica accessions, and all UK tetraploids (Cochlearia officinalis) samples/accessions for subsequent analysis in R. The subsetted UK tetraploid + UK hexaploid-only VCF produced by GATK SelectVariants is subsequently LD pruned using the prune_ld.c script written by Tuomas Hämälä (2023). 
####################

##source bash profile
source $HOME/.bash_profile

##module load gatk
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

###########
##use GATK to select all UK hexaploids and all UK tetraploids from reheadered.F4_133.ann.vcf.gz
gatk SelectVariants \
 -V ~/reheadered.F4_133.ann.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 -sn BRE_1 -sn CUM_1 -sn DAR_1 -sn DAR_3 -sn JON_001 -sn RYE_1 -sn SCO_1 -sn FRE_013 -sn FOR_1 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 \
 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 \
 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 \
 -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn GEO_2 -sn GEO_6 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL-008 \
 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_1 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 \
 -O ~/160524_data/160524_allUKtets_allUKhex.vcf.gz

##module unload
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##########

############
cd ~/160524_data/

module load samtools-uoneasy/1.18-GCC-12.3.0

##make a copy of the vcf and then unzip 
cp ./160524_allUKtets_allUKhex.vcf.gz ./160524_allUKtets_allUKhex_copy.vcf.gz

gunzip ./160524_allUKtets_allUKhex_copy.vcf.gz

module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

##########
##now ld prune the VCF file using prune_ld.c
module load gcc-uoneasy/13.2.0
##configure prune_ld
gcc ~/scripts/prune_ld.c -o ~/160524_data/160524_prune_ld -lm

#########
cd /gpfs01/home/pmyla1/160524_data/
##execute 140524_prune_ld on the 140524_UKhexaploids_only.vcf.gz
./160524_prune_ld -vcf ./160524_allUKtets_allUKhex_copy.vcf -mis 0.8 -maf 0.05 -r2 100 50 0.1 > ./160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex_copy.vcf 

module unload gcc-uoneasy/13.2.0
##########

#########
##bgzip the 140524_ld_pruned_UKhexaploids_only_copy.vcf 
cp ./160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex_copy.vcf ./160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex.vcf

module load htslib-uoneasy/1.18-GCC-13.2.0
##bgzip the ld_pruned.vcf file
bgzip ./160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex.vcf

module unload htslib-uoneasy/1.18-GCC-13.2.0
##########

echo "DONE!!"

##lastline
