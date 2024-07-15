#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=24:00:00
#SBATCH --job-name=GenotypeNOFLEET
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

mkdir -p 090624_combined_genotyped/
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
INDIR=~/300524_HaplotypeCaller_output/090624_Combined_VCF
OUTDIR=~/300524_HaplotypeCaller_output/090624_combined_genotyped/
################
##GATK GenotypeGVCFs of the additional danica samples and ionopsidium samples
 gatk GenotypeGVCFs \
   -R $REF \
   -V $INDIR/090624_Ionops_danica_NO_FLEET2.g.vcf.gz \
   -O $OUTDIR/090624_Ion.dan.NOFLEET2.genotyped.g.vcf.gz \
   -G StandardAnnotation \
   --include-non-variant-sites True
################

##unload the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
