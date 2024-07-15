#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=08:00:00
#SBATCH --job-name=Cat_R1_R2_fastq.gz
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

##cd to 170524_cutadapt/
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/

#mkdir ~/200524_new_alignments/
##use bwa mem
#bwa mem -t 20 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.sam
##now for HAM_1
#bwa mem -t 20 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_HAM_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for LWS_1
#bwa mem -t 16 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_LWS_EKDL240001890-1A_222TKYLT4_paired.sam
##now for NOT_5 
#bwa mem -t 16 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_NOT_EKDL240001890-1A_222TKYLT4_paired.sam
##now for PAR_2
#bwa mem -t 16 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_PAR_2_EKDL240001890-1A_222TKYLT4_paired.sam
##now for PEN_1
bwa mem -t 8 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_Pen_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for SBA_1
bwa mem -t 8 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_SBAY_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for SPE_2
bwa mem -t 8 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz  > ~/200524_new_alignments/200524_SPEY_2_EKDL240001890-1A_222TKYLT4_paired.sam
############

#cd ~/2024.Cochlearia.Illumina.cohort/Ionopsidium_cutadapt/
##now for the Ionopsidium accessions Iac and Ime
##first for Iac
#bwa mem -t 16 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_Iac_1P.fastq.gz ./170524_Iac_2P.fastq.gz > ~/200524_new_alignments/200524_Iac_paired.sam
##for Ime
#bwa mem -t 16 ~/C_excelsa_V5_reference/C_excelsa_V5.fa ./170524_Ime_1P.fastq.gz ./170524_Ime_2P.fastq.gz > ~/200524_new_alignments/200524_Ime_paired.sam
##unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0
#########

echo "DONE!!!"

##lastline
