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
cd /gpfs01/home/pmyla1/C_excelsa_V5_reference/

#######
##load bwa module
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

########
##use BWA to index the C_excelsa_V5.fasta
bwa index ./C_excelsa_V5.fa

##use bwa mem to align the short reads from the FRE_2 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_FRE_2.sam
##use bwa mem to align the short reads from the HAM_1 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_HAM_1.sam
####use bwa mem to align the short reads from the LWS_1 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_LWS_1.sam
##use bwa mem to align the short reads from the NOT_5 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_NOT_5.sam
##use bwa mem to align the short reads from the PAR_2 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_PAR_2.sam
##use bwa mem to align the short reads from the PEN_1 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_PEN_1.sam
##use bwa mem to align the short reads from the SBA_1 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_SBA_1.sam
##use bwa mem to align the short reads from the SPE_2 population to the C_excelsa_V5.fa 
bwa mem -t 4 ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_SPE_2.sam
#########

##unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

##########
##load samtools module for converting the sam files into bam format (-b flag) and 4 threads
module load samtools-uoneasy/1.18-GCC-12.3.0

##use samtools view to convert the alignment sam files into bam format
##first for FRE_2
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_FRE_2.sam > ../160524_new_danica_alignments/160524_aln-se_FRE_2_paired.bam
##now for HAM_1
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_HAM_1.sam > ../160524_new_danica_alignments/160524_aln-se_HAM_1_paired.bam
##next for LWS_1
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_LWS_1.sam > ../160524_new_danica_alignments/160524_aln-se_LWS_1_paired.bam
##now for NOT_5
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_NOT_5.sam > ../160524_new_danica_alignments/160524_aln-se_NOT_5_paired.bam
##now for PAR_2
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_PAR_2.sam > ../160524_new_danica_alignments/160524_aln-se_PAR_2_paired.bam
##now for PEN_1
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_PEN_1.sam > ../160524_new_danica_alignments/160524_aln-se_PEN_1_paired.bam
##now for SBA_1
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_SBA_1.sam > ../160524_new_danica_alignments/160524_aln-se_SBA_1_paired.bam
##now for SPE_2
samtools view -@ 4 -b ../160524_new_danica_alignments/160524_aln-se_SPE_2.sam > ../160524_new_danica_alignments/160524_aln-se_SPE_2_paired.bam

##########
##now use samtools sort to sort the bam files, making downstream analysis simpler
samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_FRE_2_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_FRE_2_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_HAM_1_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_HAM_1_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_LWS_1_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_LWS_1_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_NOT_5_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_NOT_5_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_PAR_2_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_PAR_2_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_PEN_1_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_PEN_1_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_SBA_1_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_SBA_1_paired.bam

samtools sort -@ 4 -o ../160524_new_danica_alignments/160524_aln-se_SPE_2_paired.sorted.bam ../160524_new_danica_alignments/160524_aln-se_SPE_2_paired.bam
###########

##########
##now use samtools index to index the final alignment files
samtools index ../160524_new_danica_alignments/160524_aln-se_FRE_2_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_HAM_1_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_LWS_1_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_NOT_5_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_PAR_2_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_PEN_1_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_SBA_1_paired.sorted.bam 

samtools index ../160524_new_danica_alignments/160524_aln-se_SPE_2_paired.sorted.bam 
##########

##########
##finally use samtools flagstat to get a summary of the overall alignment
samtools flagstat ../160524_new_danica_alignments/160524_aln-se_FRE_2_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_HAM_1_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_LWS_1_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_NOT_5_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_PAR_2_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_PEN_1_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_SBA_1_paired.sorted.bam 

samtools flagstat ../160524_new_danica_alignments/160524_aln-se_SPE_2_paired.sorted.bam 
##########


##unload samtools module
module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

echo "DONE!!!"

##lastline
