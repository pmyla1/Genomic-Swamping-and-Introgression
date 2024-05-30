#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=Cat_R1_R2_fastq.gz
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

############################
##This script was written by Luke Archer (2024) and can be used to concatenate the R1 and R2 Illumina paired-end reads (.fastq.gz) into a single "interleaved" fastq.gz file. This was performed to increase the speed of alignment production with bwa mem, because the alignments using the paired-end reads separately was extremely slow and memory intensive.
#####################

##################
##setup 
source $HOME/.bash_profile

##cd to 2024.Cochlearia.Illumina.cohort/170524_cutadapt/
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt

mkdir 180524_merged_reads/

###concatenate the R1 and R2 reads for each population into one file for alignment
cat ./170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_FLEET_2_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_180524_HAM_1_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_LWS_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_NOT_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_PAR_2_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_SBAY_1_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_SPEY_2_EKDL240001890-1A_222TKYLT4.fq.gz

cat ./170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./180524_merged_Pen_1_EKDL240001890-1A_222TKYLT4.fq.gz
###############

echo "DONE!!"

##lastline
