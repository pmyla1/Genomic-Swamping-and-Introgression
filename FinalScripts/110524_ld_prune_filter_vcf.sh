#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=28g
#SBATCH --time=02:00:00
#SBATCH --job-name=LD_prune_filter_vcf
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load samtools module to unzip the vcf file
#module load samtools-uoneasy/1.18-GCC-12.3.0

##unzip the VCF file you want to LD prune
#gunzip /gpfs01/home/pmyla1/bcftools_gatk_output/110524_WG_allUKhex_allUKdips_someUKtets.vcf.gz

##unload samtools 
#module unload samtools-uoneasy/1.18-GCC-12.3.0

##load gcc to compile prune_ld.c script by Tuomas Hamala (2024)
#module load gcc-uoneasy/13.2.0

###### 
##compile prune_ld.c script
#gcc /gpfs01/home/pmyla1/scripts/prune_ld.c -o /gpfs01/home/pmyla1/bcftools_gatk_output/prune_ld -lm 

##execute the script on the whole genome UK_dips_tets_danica_anglica.vcf.gz
#/gpfs01/home/pmyla1/bcftools_gatk_output/prune_ld -vcf /gpfs01/home/pmyla1/bcftools_gatk_output/110524_WG_allUKhex_allUKdips_someUKtets.vcf -mis 0.9 -maf 0.05 -r2 100 50 0.1 > /gpfs01/home/pmyla1/bcftools_gatk_output/ld_pruned_110524_WG_allUKhex_allUKdips_someUKtets.vcf
 
##unload gcc module 
#module unload gcc-uoneasy/13.2.0

##load htslib
module load htslib-uoneasy/1.18-GCC-13.2.0

#bgzip the newly produced ld_pruned vcf 
#bgzip /gpfs01/home/pmyla1/bcftools_gatk_output/ld_pruned_110524_WG_allUKhex_allUKdips_someUKtets.vcf

##bgzip the 110524_WG_allUKhex_allUKdips_someUKtets.vcf
bgzip /gpfs01/home/pmyla1/bcftools_gatk_output/110524_WG_allUKhex_allUKdips_someUKtets.vcf

##unload htslib
module unload htslib-uoneasy/1.18-GCC-13.2.0
######

echo 'DONE!!'

##lastline
