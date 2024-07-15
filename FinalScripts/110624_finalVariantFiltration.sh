#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=24:00:00
#SBATCH --job-name=SelectVariantsDepthFilter
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/

###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F2.best.practice.g.vcf.gz
OUTMASK=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_depth.mask.Ion.dan.g.vcf.gz
OUTF3=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_Ion.dan.filtered.F3.g.vcf.gz
OUTF4=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_Ion.dan.filtered.F4.g.vcf.gz
################
##GATK SelectVariants to select biallelic variants only
gatk SelectVariants \
   -R $REF \
   -V $VCF \
   -O $OUTMASK \
   --select "DP<136" ##depth cut off = 1.6*mean depth
###############
gatk VariantFiltration \
    -R $REF \
    -V $VCF \
    -O $OUTF3 \
    --mask $OUTMASK \
    --filter-not-in-mask
###############
gatk SelectVariants \
    -R $REF \
    -V $OUTF3 \
    -O $OUTF4 \
    --exclude-filtered True
#############

##unload the GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
