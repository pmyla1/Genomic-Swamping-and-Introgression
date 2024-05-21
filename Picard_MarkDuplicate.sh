#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=Picard_mark_duplicates
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

###############################
##This script was written by Luke Archer (2024) and uses Picard MarkDuplicates (version 3.3.0) to mark the duplicate reads in the .sorted.bam files produced by convert_to_bam_samtools_stats.sh, and produces duplicate-marked .sorted.bam files. 
##############################

##################
##setup 
source $HOME/.bash_profile

############
##load the picard module to mark the duplicates in the bam alignment files
module load picard-uoneasy/3.0.0-Java-17

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

##use MarkDuplicates to mark the duplicated reads in FLE_2 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_FLE_2_paired.sorted.bam \
      O=./190524_FLE_2_marked_duplicates.sorted.bam \
      M=./190524_FLE_2_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in HAM_1 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_HAM_1_paired.sorted.bam \
      O=./190524_HAM_1_marked_duplicates.sorted.bam \
      M=./190524_HAM_1_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in LWS_1 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_LWS_1_paired.sorted.bam \
      O=./190524_LWS_1_marked_duplicates.sorted.bam \
      M=./190524_LWS_1_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in NOT_5 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_NOT_5_paired.sorted.bam \
      O=./190524_NOT_5_marked_duplicates.sorted.bam \
      M=./190524_NOT_5_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in PAR_2 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_PAR_2_paired.sorted.bam \
      O=./190524_PAR_2_marked_duplicates.sorted.bam \
      M=./190524_PAR_2_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in PEN_1 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_PEN_1_paired.sorted.bam \
      O=./190524_PEN_1_marked_duplicates.sorted.bam \
      M=./190524_PEN_1_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in SBA_1 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_SBA_1_paired.sorted.bam \
      O=./190524_SBA_1_marked_duplicates.sorted.bam \
      M=./190524_SBA_1_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in SPE_2 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_SPE_2_paired.sorted.bam \
      O=./190524_SPE_2_marked_duplicates.sorted.bam \
      M=./190524_SPE_2_marked_duplicates_metrics.txt
#########

############
##now do the same for the Ionopsidium samples
cd
~/2024.Cochlearia.Illumina.cohort/Ionopsidium_cutadapt/180524_Ionopsidium_merged_reads/180524_Ionopsidium_alignments/190524_bam_files/
############
##use MarkDuplicate to mark the duplicated reads in SBA_1 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_Iac_paired.sorted.bam \
      O=./190524_Iac_marked_duplicates.sorted.bam \
      M=./180524_Iac_marked_duplicates_metrics.txt
#########
##use MarkDuplicate to mark the duplicated reads in SPE_2 
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=./180524_Ime_paired.sorted.bam \
      O=./190524_Ime_marked_duplicates.sorted.bam \
      M=./190524_Ime_marked_duplicates_metrics.txt
#########

###unload module
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline
