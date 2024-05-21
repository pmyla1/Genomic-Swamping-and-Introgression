#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --job-name=CollectSummaryAlignmentMetrics
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the java picard module to grep the read groups from the duplicate marked BAM files
module load picard-uoneasy/3.0.0-Java-17 

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

mkdir ~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/
#############
##calculate alignment summary metrics for FLE_2 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_FLE_2_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_FLE_2_aln_sum.txt
############
##calculate alignment summary metrics for HAM_1 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_HAM_1_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_HAM_1_aln_sum.txt
############
##calculate alignment summary metrics for LWS_1 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_LWS_1_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_LWS_1_aln_sum.txt
############
##calculate alignment summary metrics for NOT_5 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_NOT_5_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_NOT_5_aln_sum.txt
############
##calculate alignment summary metrics for PAR_2 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_PAR_2_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_PAR_2_aln_sum.txt
############
##calculate alignment summary metrics for PEN_1 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_PEN_1_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_PEN_1_aln_sum.txt
#############
##calculate alignment summary metrics for SBA_1 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_SBA_1_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_SBA_1_aln_sum.txt
#############
##calculate alignment summary metrics for SPE_2 duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_SPE_2_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_SPE_2_aln_sum.txt
#############
##calculate alignment summary metrics for Iac duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_Iac_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_Iac_aln_sum.txt
##############
##calculate alignment summary metrics for Ime duplicate marked bam
    java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics \
          R=~/C_excelsa_V5_reference/C_excelsa_V5.fa \
          I=./190524_Ime_marked_duplicates.sorted.bam \
          O=~/2024.Cochlearia.Illumina.cohort/200524_Alignment_Summary_Metrics/200524_Ime_aln_sum.txt
##############

module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline

