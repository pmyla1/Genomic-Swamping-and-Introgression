#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=ValidateSamFile
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

#############################
##This script was written by Luke Archer (2024) and uses Picard ValidateSamFile (version 3.3.0) to validate the duplicate marked .sorted.bam files produced by Picard_MarkDuplicate.sh. The output is a "validated.txt" file which indicates that the .sorted.bam files do not have a read group associated with them, therefore GATK HaplotypeCaller cannot call haplotypes on these files.
############################

##################
##setup 
source $HOME/.bash_profile

############
##load the picard module 
module load picard-uoneasy/3.0.0-Java-17

##change directory to the 190524_bam_alignments
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/
############
##validate each sam file
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_FLE_2_marked_duplicates.sorted.bam \
      -O ./200524_FLE_2_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_HAM_1_marked_duplicates.sorted.bam \
      -O ./200524_HAM_1_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_LWS_1_marked_duplicates.sorted.bam \
      -O ./200524_LWS_1_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_NOT_5_marked_duplicates.sorted.bam \
      -O ./200524_NOT_5_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_PAR_2_marked_duplicates.sorted.bam \
      -O ./200524_PAR_2_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_PEN_1_marked_duplicates.sorted.bam \
      -O ./200524_PEN_1_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_SBA_1_marked_duplicates.sorted.bam \
      -O ./200524_SBA_1_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_SPE_2_marked_duplicates.sorted.bam \
      -O ./200524_SPE_2_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_Iac_marked_duplicates.sorted.bam \
      -O ./200524_Iac_validated.txt \
      --MODE SUMMARY
############
java -jar $EBROOTPICARD/picard.jar ValidateSamFile \
      -I ./190524_Ime_marked_duplicates.sorted.bam \
      -O ./200524_Ime_validated.txt \
      --MODE SUMMARY
############

###unload module
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline
