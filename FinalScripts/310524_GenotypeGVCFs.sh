#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40g
#SBATCH --time=08:00:00
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

mkdir ~/310524_Genotyped_VCF/

###make environmental variables for the reference genome and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
INVCF=~/300524_HaplotypeCaller_output/Combined_Ionopsidium_VCF/300524_combined_Ionopsidium.g.vcf.gz
OUTDIR=~/310524_Genotyped_VCF
################
##GATK GenotypeGVCFs 
 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $REF \
   -V $INVCF \
   -O $OUTDIR/310524_Genotyped_Ionopsidium.vcf.gz \
   -G StandardAnnotation \
   --include-non-variant-sites true \
   --sample-ploidy 32
##############   

##unload the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
