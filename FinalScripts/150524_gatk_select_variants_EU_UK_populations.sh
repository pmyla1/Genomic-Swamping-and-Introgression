#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=16g
#SBATCH --time=01:00:00
#SBATCH --job-name=gatk_select_from_all_areas
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##module load gatk
#module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

###########
##use GATK to select only the UK hexaploids from the original VCF file
#gatk SelectVariants \
# -V /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz \
# --select-type-to-include SNP \
# --restrict-alleles-to BIALLELIC -xl-sn ROT_004 -xl-sn ROT_006 -xl-sn ROT_007 -xl-sn ROT_013 -xl-sn ELI_001 -xl-sn ELI_002 -xl-sn ELI_003 -xl-sn ELI_004 -xl-sn LNL_001 -xl-sn LNL_002 -xl-sn LNL_003 -xl-sn LNL_008 -xl-sn NEI_8 -xl-sn NEI_9 -xl-sn NEI_1 -xl-sn NEI_2 -xl-sn NEI_3 -xl-sn NEI_4 -xl-sn GEO_1 -xl-sn GEO_2 -xl-sn GEO_6 -xl-sn BNK_21 -xl-sn CHA_1 -xl-sn CHA_2 -xl-sn SKI_004 -xl-sn SKI_005 -xl-sn MEL_001 -xl-sn MEL_002 -xl-sn BEA_002 -xl-sn BEA_004 -xl-sn BEA_010 -xl-sn RUZ_001 -O /gpfs01/home/pmyla1/150524_VCFs_data/150524_EU_UK_pops.vcf.gz

##module unload
#module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##########

############
#cd /gpfs01/home/pmyla1/150524_VCFs_data/

#module load samtools-uoneasy/1.18-GCC-12.3.0
##make a copy of the vcf and then unzip 
#cp ./150524_EU_UK_pops.vcf.gz ./150524_EU_UK_pops_copy.vcf.gz

#gunzip ./150524_EU_UK_pops_copy.vcf.gz

#module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

##########
##now ld prune the VCF file using prune_ld.c
#module load gcc-uoneasy/13.2.0
##configure prune_ld
#gcc /gpfs01/home/pmyla1/scripts/prune_ld.c -o /gpfs01/home/pmyla1/150524_VCFs_data/150524_prune_ld -lm

#########
#cd /gpfs01/home/pmyla1/150524_VCFs_data/
##execute 140524_prune_ld on the 140524_UKhexaploids_only.vcf.gz
#./150524_prune_ld -vcf ./150524_EU_UK_pops_copy.vcf -mis 0.8 -maf 0.05 -r2 100 50 0.1 > ./150524_ld_pruned_20PCTmis_maf005_EU_UK_pops_copy.vcf 

#module unload gcc-uoneasy/13.2.0
##########

#########
##bgzip the 140524_ld_pruned_UKhexaploids_only_copy.vcf 
#cp ./150524_ld_pruned_20PCTmis_maf005_EU_UK_pops_copy.vcf ./150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf

module load htslib-uoneasy/1.18-GCC-13.2.0
##bgzip the ld_pruned.vcf file
bgzip /gpfs01/home/pmyla1/150524_VCFs_data/150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf

module unload htslib-uoneasy/1.18-GCC-13.2.0
##########

echo "DONE!!"

##lastline
