#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=BWA_align_reads
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##cd to pmyla1/
cd /gpfs01/home/pmyla1/

#######
##load bwa module
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

########
##use BWA to index the C_excelsa_V5.fasta, using the bwtsw algorithim used for the human genome
bwa index -p 160524_C_excelsa_V5 -a bwtsw ./C_excelsa_V5.fa

##use bwa mem to align the short reads from the FRE_2 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_FRE_2.sam
##use bwa mem to align the short reads from the HAM_1 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_HAM_1.sam
####use bwa mem to align the short reads from the LWS_1 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/LWS_1/LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/LWS_1/LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_LWS_1.sam
##use bwa mem to align the short reads from the NOT_5 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/NOT_5/NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/NOT_5/NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_NOT_5.sam
##use bwa mem to align the short reads from the PAR_2 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/PAR_2/PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/PAR_2/PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_PAR_2.sam
##use bwa mem to align the short reads from the PEN_1 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/PEN_1/Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/PEN_1/Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_PEN_1.sam
##use bwa mem to align the short reads from the SBA_1 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/SBA_1/SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/SBA_1/SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_SBA_1.sam
##use bwa mem to align the short reads from the SPE_2 population to the C_excelsa_V5.fa 
bwa mem -P 160524_C_excelsa_V5 ./C_excelsa_V5.fa ./2024.Cochlearia.Illumina.cohort/SPE_2/SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/SPE_2/SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ./160524_new_danica_alignments/160524_aln-se_SPE_2.sam
#########

##unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

echo "DONE!!!"

##lastline