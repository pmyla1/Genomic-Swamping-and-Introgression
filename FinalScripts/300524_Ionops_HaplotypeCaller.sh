#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40g
#SBATCH --time=10:00:00
#SBATCH --job-name=Ionops_HaplotypeCaller
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the GATK module for haplotype calling, specifying -ploidy 6 
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

########
##make a sequence dictionary for the C_excelsa_V5 reference genome
gatk CreateSequenceDictionary -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa

module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
########

#########
##load samtools to index fasta file
module load samtools-uoneasy/1.18-GCC-12.3.0

samtools faidx ~/C_excelsa_V5_reference/C_excelsa_V5.fa

module unload samtools-uoneasy/1.18-GCC-12.3.0
#########

#############
##re-load gatk module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##change directory to the duplicate marked bams
cd ~/ionops/04.marked/

mkdir ~/300524_HaplotypeCaller_output/

##make environmental variables for the reference genome and for the output directory
OUTDIR=~/300524_HaplotypeCaller_output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa

#########
##use GATK HaplotypeCaller to make a gVCF of the additional outgroup Ionopsidium files 
##firstly for Iac
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Iac/Iac.mrkd.nmd.bam \
#   -O $OUTDIR/Iac.g.vcf.gz \
#   -bamout $OUTDIR/Iac.raw.8x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 8 \
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
##now for Ime
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./Ime/Ime.mrkd.nmd.bam \
   -O $OUTDIR/Ime.g.vcf.gz \
   -bamout $OUTDIR/Ime.raw.8x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 8 \
############
##now for Iab_1
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Iab_1/Iac_1.mrkd.nmd.bam \
#   -O $OUTDIR/Iab_1.g.vcf.gz \
#   -bamout $OUTDIR/Iab_1.raw.8x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 8 \
############
##now for Iab_2
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Iab_2/Iac_2.mrkd.nmd.bam \
#   -O $OUTDIR/Iab_2.g.vcf.gz \
#   -bamout $OUTDIR/Iab_2.raw.8x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 8 \
############

###unload module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
