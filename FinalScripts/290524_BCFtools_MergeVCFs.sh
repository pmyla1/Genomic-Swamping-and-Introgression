#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --job-name=BCFtoolsMergeVCFs
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##module load
module load bcftools-uoneasy/1.18-GCC-13.2.0

cd /gpfs01/home/pmyla1/

##############
##use bcftools to merge reheadered.F4_133.ann.vcf.gz with the filtered.F4.vcf.gz
bcftools merge --print-header -Oz -o ./290524_Ionops_Cochlearia_merged.vcf.gz ./reheadered.F4_133.ann.vcf.gz ./filtered.F4.vcf.gz


##unload module
module unload bcftools-uoneasy/1.18-GCC-13.2.0
