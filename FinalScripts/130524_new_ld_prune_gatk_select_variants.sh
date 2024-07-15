#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=28g
#SBATCH --time=02:00:00
#SBATCH --job-name=sele_all_danica_anglica_UK_dips
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

#####
##index the reheadered.F4_133.ann.vcf.gz file with gatk IndexFeatureFile
gatk IndexFeatureFile -I ~/reheadered.F4_133.ann.vcf.gz

######
##use GATK SelectVariants to select all UK diploids, tetraploids, C. danica, and putative C. anglica
gatk SelectVariants -V ~/reheadered.F4_133.ann.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 \
 -sn BNK_21 -sn BRE_1 -sn CHA_1 -sn CHA_2 -sn CUM_1 \
 -sn DAR_1 -sn DAR_3 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 \
 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FOR_1 -sn FRE_013 \
 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn GEO_2 -sn GEO_6 \
 -sn JON_001 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 \
 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 \
 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 \
 -sn LOS_1 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 \
 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 \
 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 \
 -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn RYE_1 -sn SCO_1 \
 -sn SCU_1 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 \
 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 \
 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 -O ~/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex.vcf.gz \
 --allow-nonoverlapping-command-line-samples
##########

##unload GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

############
##now use samtools to unzip the newly produced VCF 
##load samtools module to unzip the vcf file
module load samtools-uoneasy/1.18-GCC-12.3.0

##make a copy of the 130524_allUKdips_allUKtets_allUKhex.vcf.gz
cp ~/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex.vcf.gz /gpfs01/home/pmyla1/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex_copy.vcf.gz

##unzip the VCF file you want to LD prune
gunzip ~/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex_copy.vcf.gz

##unload samtools 
module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

##########
##Compile and execute prune_ld.c script
##load gcc to compile prune_ld.c script by Tuomas Hamala (2024)
module load gcc-uoneasy/13.2.0

###### 
##compile prune_ld.c script
gcc ~/scripts/prune_ld.c -o ~/bcftools_gatk_output/prune_ld -lm 

##execute the script on the whole genome UK_dips_tets_danica_anglica.vcf
~/bcftools_gatk_output/prune_ld -vcf ~/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex_copy.vcf -mis 0.9 -maf 0.05 -r2 100 50 0.1 > /gpfs01/home/pmyla1/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf
 
##unload gcc module 
module unload gcc-uoneasy/13.2.0

##load htslib
module load htslib-uoneasy/1.18-GCC-13.2.0

#bgzip the newly produced ld_pruned vcf 
bgzip ~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf

##bgzip the 110524_WG_allUKhex_allUKdips_someUKtets.vcf
bgzip ~/bcftools_gatk_output/130524_allUKdips_allUKtets_allUKhex.vcf

##unload htslib
module unload htslib-uoneasy/1.18-GCC-13.2.0

echo "DONE!!"

##lastline
