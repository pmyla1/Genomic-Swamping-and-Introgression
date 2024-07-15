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
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

#####
##index the reheadered.F4_133.ann.vcf.gz file with gatk IndexFeatureFile
gatk IndexFeatureFile -I /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz

######
##use GATK SelectVariants to select all UK diploids, tetraploids, C. danica, and putative C. anglica
gatk SelectVariants -V /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --select "AF > 0.5 && AF < 1.0" \
 -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 \
 -sn BNK_21 -sn BRE_1 -sn CHA_1 -sn CHA_2 -sn CUM_1 \
 -sn DAR_1 -sn DAR_3 \
 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FOR_1 -sn FRE_013 \
 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 \
 -sn JON_001 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 \
 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 \
 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 \
 -sn LOS_1 -sn LOS_6 -sn LOS_7 \
 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 \
 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 \
 -sn RYE_1 -sn SCO_1 \
 -sn SCU_1 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 \
 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 \
 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 -O /gpfs01/home/pmyla1/bcftools_gatk_output/110524_WG_allUKhex_allUKdips_someUKtets.vcf.gz \
 --allow-nonoverlapping-command-line-samples
##########

##unload GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
