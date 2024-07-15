#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=08:00:00
#SBATCH --job-name=SBAY-SPEY
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile
##load BWA for alignments
#module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

#make a variable for the "metadata" 
metadata=EKDL240001890-1A_222TKYLT4

#store the output directory and reference genome in variables called OUT & REF, respectively
OUT=~/220524_alignments
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
####################
##cd to the fastq.gz folder
#cd ~/2024.Cochlearia.Illumina.cohort/220524_trimmed_reads/

##############
#bwa mem \
#     -t 16 $REF \
#     ./SBAY_1_${metadata}_L1_1.trimmed.fq.gz ./SBAY_1_${metadata}_L1_2.trimmed.fq.gz \
#     > $OUT/SBAY_1_${metadata}_aln-pe.sam
##############
#bwa mem \
#     -t 16 $REF \
#     ./SPEY_2_${metadata}_L1_1.trimmed.fq.gz ./SPEY_2_${metadata}_L1_2.trimmed.fq.gz \
#     > $OUT/SPEY_2_${metadata}_aln-pe.sam

#module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0
############
##Samtools to sort, index, and alignment summary
module load samtools-uoneasy/1.18-GCC-12.3.0

cd ~/220524_alignments/
#on SBA_1
samtools view -@ 4 -h -b ./SBAY_1_${metadata}_aln-pe.sam -o ./bam_files/SBAY_1_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/SBAY_1_${metadata}.sorted.bam ./bam_files/SBAY_1_${metadata}.bam
samtools index ./bam_files/SBAY_1_${metadata}.sorted.bam
samtools flagstat ./bam_files/SBAY_1_${metadata}.sorted.bam > ./bam_files/SBAY_1_${metadata}.flagstats
############
##SPE_2
samtools view -@ 4 -h -b ./SPEY_2_${metadata}_aln-pe.sam -o ./bam_files/SPEY_2_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/SPEY_2_${metadata}.sorted.bam ./bam_files/SPEY_2_${metadata}.bam 
samtools index ./bam_files/SPEY_2_${metadata}.sorted.bam
samtools flagstat ./bam_files/SPEY_2_${metadata}.sorted.bam > ./bam_files/SPEY_2_${metadata}.flagstats
##########

module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
