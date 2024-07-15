#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8g
#SBATCH --time=01:00:00
#SBATCH --job-name=MARKDUPLICATES
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the Picard module
module load picard-uoneasy/3.0.0-Java-17

cd ~/220524_alignments/bam_files/

mkdir -p duplicate_marked_bams/

##make environmental variables
OUTDIR=~/220524_alignments/bam_files/duplicate_marked_bams
meta=EKDL240001890-1A_222TKYLT4
##first on FLEET_2
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./FLEET_2_${meta}.sorted.bam -O $OUTDIR/FLEET_2_${meta}.marked_duplicates.bam -M $OUTDIR/FLEET_2_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
############
##try on Iac
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./Iac.sorted.bam -O $OUTDIR/Iac.marked_duplicates.bam -M $OUTDIR/Iac.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
#########
##try on Ime
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./Ime.sorted.bam -O $OUTDIR/Ime.marked_duplicates.bam -M $OUTDIR/Ime.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
#########
##now on LWS
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./LWS_${meta}.sorted.bam -O $OUTDIR/LWS_${meta}.marked_duplicates.bam -M $OUTDIR/LWS_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
#########
##now on NOT
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./NOT_${meta}.sorted.bam -O $OUTDIR/NOT_${meta}.marked_duplicates.bam -M $OUTDIR/NOT_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
#########
##now on PAR_2
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./PAR_2_${meta}.sorted.bam -O $OUTDIR/PAR_2_${meta}.marked_duplicates.bam -M $OUTDIR/PAR_2_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
##########
##now on SBAY_1
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./SBAY_1_${meta}.sorted.bam -O $OUTDIR/SBAY_1_${meta}.marked_duplicates.bam -M $OUTDIR/SBAY_1_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
###########
##now on SPEY_2
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./SPEY_2_${meta}.sorted.bam -O $OUTDIR/SPEY_2_${meta}.marked_duplicates.bam -M $OUTDIR/SPEY_2_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
###########
##now on HAM_1
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./HAM_1_${meta}.sorted.bam -O $OUTDIR/HAM_1_${meta}.marked_duplicates.bam -M $OUTDIR/HAM_1_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
############
##now on Pen_1
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./Pen_1_${meta}.sorted.bam -O $OUTDIR/Pen_1_${meta}.marked_duplicates.bam -M $OUTDIR/Pen_1_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true
###########

module load picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline
