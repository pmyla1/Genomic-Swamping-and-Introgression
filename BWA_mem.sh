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
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

##cd to 170524_cutadapt/
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/

mkdir 180524_alignments/
##use bwa mem on one of the merged fastq.gz files to align to the reference
##first for FLE_2
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_FLEET_2_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.sam
##now for HAM_1
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_180524_HAM_1_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_180524_HAM_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for LWS_1
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_LWS_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_LWS_EKDL240001890-1A_222TKYLT4_paired.sam
##now for NOT_5 
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_NOT_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_NOT_EKDL240001890-1A_222TKYLT4_paired.sam
##now for PAR_2
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_PAR_2_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_PAR_2_EKDL240001890-1A_222TKYLT4_paired.sam
##now for PEN_1
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_Pen_1_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_Pen_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for SBA_1
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_SBAY_1_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_SBAY_1_EKDL240001890-1A_222TKYLT4_paired.sam
##now for SPE_2
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_SPEY_2_EKDL240001890-1A_222TKYLT4.fq.gz  > ./180524_alignments/180524_SPEY_2_EKDL240001890-1A_222TKYLT4_paired.sam
############

##unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

echo "DONE!!"

##lastline


