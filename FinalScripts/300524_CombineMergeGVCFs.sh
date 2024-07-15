#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=CombineGVCFs
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

cd ~/300524_HaplotypeCaller_output/

mkdir -p Combined_Ionopsidium_VCF/
###make environmental variables for the reference genome and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
OUTDIR=~/300524_HaplotypeCaller_output/Combined_Ionopsidium_VCF
################
##GATK CombineGVCFs 
 gatk CombineGVCFs \
   -R $REF \
   --variant ./Iab_1.g.vcf.gz \
   --variant ./Iab_2.g.vcf.gz \
   --variant ./Iac.g.vcf.gz \
   --variant ./Ime.g.vcf.gz \
   -O $OUTDIR/300524_combined_Ionopsidium.g.vcf.gz
################

##unload the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
