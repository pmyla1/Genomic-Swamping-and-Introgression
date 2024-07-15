#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=24g
#SBATCH --time=01:00:00
#SBATCH --job-name=sele_fewer_individuals
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load GATK module
#module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

#####
##index the reheadered.F4_133.ann.vcf.gz file with gatk IndexFeatureFile
#gatk IndexFeatureFile -I ~/reheadered.F4_133.ann.vcf.gz

######
##use GATK SelectVariants to select all UK diploids, tetraploids, C. danica, and putative C. anglica
#gatk SelectVariants -V ~/reheadered.F4_133.ann.vcf.gz \
# --select-type-to-include SNP \
# --restrict-alleles-to BIALLELIC \
# -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
# -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 \
# -sn BRE_1 -sn CUM_1 \
# -sn DAR_1 -sn DAR_3 \
# -sn FOR_1 -sn FRE_013 \
# -sn JON_001 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 \
# -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 \
# -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 \
# -sn GEO_2 -sn GEO_6 \
# -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 \
# -sn RYE_1 -sn SCO_1 \
# -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 \
# -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
# -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
# -O ~/050624_VCF/050624_WG_allUKhex_someUKdips_someUKtets.vcf.gz \
# --allow-nonoverlapping-command-line-samples
##########

##unload GATK module
#module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##############
##unzip the 050624_WG_allUKhex_someUKdips_someUKtets.vcf.gz with gunzip 
##now use samtools to unzip the newly produced VCF 
##load samtools module to unzip the vcf file
#module load samtools-uoneasy/1.18-GCC-12.3.0

##unzip the VCF file you want to LD prune
#gunzip ~/050624_VCF/050624_WG_allUKhex_someUKdips_someUKtets.vcf.gz

##unload samtools 
#module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

##########
##Compile and execute prune_ld.c script
##load gcc to compile prune_ld.c script by Tuomas Hamala (2024)
#module load gcc-uoneasy/13.2.0

#gcc ~/scripts/prune_ld.c -o ~/scripts/prune_ld -lm
###### 
##LD prune the 050624_WG_allUKhex_someUKdips_someUKtets.vcf.gz
#~/scripts/prune_ld -vcf ~/050624_VCF/050624_WG_allUKhex_someUKdips_someUKtets.vcf -mis 0.9 -maf 0.05 -r2 100 50 0.1 > ~/050624_VCF/050624_ld_pruned_allUKhex_someUKdips_someUKtets.vcf
 
##unload gcc module 
#module unload gcc-uoneasy/13.2.0

##############
##load htslib to use bzgip
module load htslib-uoneasy/1.18-GCC-13.2.0

#bgzip the newly produced ld_pruned vcf 
bgzip ~/050624_VCF/050624_ld_pruned_allUKhex_someUKdips_someUKtets.vcf

##unload htslib
module unload htslib-uoneasy/1.18-GCC-13.2.0

echo "DONE!!"

##lastline

echo "DONE!!"

##lastline
