#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20g
#SBATCH --time=04:00:00
#SBATCH --job-name=SelectALLVariants
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

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
 -O ~/050624_VCF/060624_AllSamples_BISNP.vcf.gz
#########
##module unload
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##########

##################
#SAMTOOLS TO UNZIP THE VCF FILE IN PREPARATION FOR PRUNE_LD.C 
############
cd ~/050624_VCF/

module load samtools-uoneasy/1.18-GCC-12.3.0
##make a copy of the vcf and then unzip 
cp ./060624_AllSamples_BISNP.vcf.gz ./060624_AllSamples_BISNP_copy.vcf.gz

#unzip the copy
gunzip ./060624_AllSamples_BISNP_copy.vcf.gz

module unload samtools-uoneasy/1.18-GCC-12.3.0
##################

#################
##LD PRUNE THE VCF FILE USING PRUNE_LD.C SCRIPT WRITTEN BY HÄMÄLÄ (2024)
module load gcc-uoneasy/13.2.0
##configure prune_ld
#gcc ~/scripts/prune_ld.c -o ~/scripts/prune_ld -lm

#########
##execute 140524_prune_ld on the 060624_AllSamples_BISNP_copy.vcf - maximum missing data (10%), minor allele freq (0.05), squared correlation 
~/scripts/prune_ld -vcf ./060624_AllSamples_BISNP_copy.vcf -mis 0.9 -maf 0.05 -r2 50 10 0.1 > ./060524_LDpruned_allSamples_BiallelicSNPs.vcf  

module unload gcc-uoneasy/13.2.0
##########

#########
##bgzip the 140524_ld_pruned_UKhexaploids_only_copy.vcf 
cp ./060524_LDpruned_allSamples_BiallelicSNPs.vcf ./060524_LDpruned_allSamples_BiallelicSNPs_copy.vcf

module load htslib-uoneasy/1.18-GCC-13.2.0
##bgzip the ld_pruned.vcf file
bgzip ./060524_LDpruned_allSamples_BiallelicSNPs.vcf

module unload htslib-uoneasy/1.18-GCC-13.2.0
##########

echo "DONE!!"

##lastline
