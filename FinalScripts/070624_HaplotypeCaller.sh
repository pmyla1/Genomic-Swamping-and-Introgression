#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=40g
#SBATCH --time=16:00:00
#SBATCH --job-name=HaplotypeCaller
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
cd ~/220524_alignments/bam_files/duplicate_marked_bams_RG/

mkdir ~/300524_HaplotypeCaller_output/

##make environmental variables for the reference genome and for the output directory
OUTDIR=~/300524_HaplotypeCaller_output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
###########
##for FLEET_2
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./FLEET_2_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
#   -O $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
#   -bamout $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.6x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for LWS_1
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./LWS_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
   -O $OUTDIR/LWS_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -bamout $OUTDIR/LWS_EKDL240001890-1A_222TKYLT4.6x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 6 \
############
##now for NOT
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./NOT_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
   -O $OUTDIR/NOT_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -bamout $OUTDIR/NOT_EKDL240001890-1A_222TKYLT4.6x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 6 \
############
##now for PAR_2
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./PAR_2_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
   -O $OUTDIR/PAR_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -bamout $OUTDIR/PAR_2_EKDL240001890-1A_222TKYLT4.6x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 6 \
############
##now for PEN_1
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Pen_1_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
#   -O $OUTDIR/Pen_1_EKDL240001890-1A_222TKYLT4.g.vcf.gz \#
#   -bamout $OUTDIR/Pen_1_EKDL240001890-1A_222TKYLT4.6x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for SBAY_1
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./SBAY_1_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
#   -O $OUTDIR/SBAY_1_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
#   -bamout $OUTDIR/SBAY_1_EKDL240001890-1A_222TKYLT4.6x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for SPEY_2
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./SPEY_2_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
#   -O $OUTDIR/SPEY_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
#   -bamout $OUTDIR/SPEY_2_EKDL240001890-1A_222TKYLT4.6x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for Iac
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Iac.marked_duplicates.bam \
#   -O $OUTDIR/Iac.g.vcf.gz \
#   -bamout $OUTDIR/Iac.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for Ime
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./Ime.marked_duplicates.bam \
#   -O $OUTDIR/Ime.g.vcf.gz \
#   -bamout $OUTDIR/Ime.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
##now for HAM_1
# gatk --java-options "-Xmx4g" HaplotypeCaller  \
#   -R $REF \
#   -I ./HAM_1_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
#   -O $OUTDIR/HAM_1_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
#   -bamout $OUTDIR/HAM_1_EKDL240001890-1A_222TKYLT4.6x.bam \
#   --emit-ref-confidence BP_RESOLUTION \
#   --min-base-quality-score 25 \
#   --minimum-mapping-quality 25 \
#   --sample-ploidy 6 \
############
############

###unload module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
