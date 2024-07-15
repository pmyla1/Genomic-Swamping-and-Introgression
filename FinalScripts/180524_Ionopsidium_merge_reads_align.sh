#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=Cat_align_Ionopsidium
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

###change working directory to the Ionopsidium cutadapt reads
#cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/Ionopsidium_cutadapt/

##make a new directory to store the output
#mkdir 180524_Ionopsidium_merged_reads
###concatenate the 1P & 2P reads for both Ime and Iac
#cat ./170524_Iac_1P.fastq.gz ./170524_Iac_2P.fastq.gz > ./180524_Ionopsidium_merged_reads/180524_Iac_merged.fastq.gz

#cat ./170524_Ime_1P.fastq.gz ./170524_Ime_2P.fastq.gz > ./180524_Ionopsidium_merged_reads/180524_Ime_merged.fastq.gz
###############

#############
##now use bwa mem to align the Ionopsidium reads to the C_excelsa_V5.fa reference

##firstly change directory to the 180524_Ionopsidium_merged_reads/
#cd 180524_Ionopsidium_merged_reads/

##make a new directory to store the alignments
#mkdir 180524_Ionopsidium_alignments/

##load bwa module
#module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

###########
##firstly align the Iac merged reads onto the C_excelsa_V5.fa reference
#bwa mem -t 4 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_Iac_merged.fastq.gz > ./180524_Ionopsidium_alignments/180524_Iac_paired.sam

##now do the same for Ime merged reads
#bwa mem -t 4 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_Ime_merged.fastq.gz > ./180524_Ionopsidium_alignments/180524_Ime_paired.sam
##########

##unload bwa module
#module unload module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

##########
##change directory to the 180524_Ionopsidium_alignments
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/Ionopsidium_cutadapt/180524_Ionopsidium_merged_reads/180524_Ionopsidium_alignments/
##load samtools view to convert the sam files to bam files
module load samtools-uoneasy/1.18-GCC-12.3.0

##convert Iac sam to bam 
#samtools view -@ 8 -b ./180524_Iac_paired.sam > ./180524_Iac_paired.bam
##use samtools sort to produce a sorted Iac bam file
#samtools sort -@ 8 -o ./180524_Iac_paired.sorted.bam ./180524_Iac_paired.bam
##finally index the sorted Iac bam file with samtools index
#samtools index ./180524_Iac_paired.sorted.bam
##now get a summary of the alignment with samtools flagstat
samtools flagstat ./180524_Iac_paired.sorted.bam > ./180524_Iac_alignment_summary.txt

##convert Ime sam to bam
#samtools view -@ 8 -b ./180524_Ime_paired.sam > ./180524_Ime_paired.bam
##use samtools sort to produce a sorted Ime bam file
#samtools sort -@ 8 -o ./180524_Ime_paired.sorted.bam ./180524_Ime_paired.bam 
#finally index the sorted Ime bam file with samtools index
#samtools index ./180524_Ime_paired.sorted.bam
##now get a summary of the alignment with samtools flagstat
samtools flagstat ./180524_Ime_paired.sorted.bam > ./180524_Ime_alignment_summary.txt

module unload samtools-uoneasy/1.18-GCC-12.3.0
###########

#############
##make a new directory to store the different output file types
mkdir 190524_alignment_summaries

mkdir 190524_bam_files

mkdir 190524_sam_files
###########

############
##move the different file types to the appropriate folder
mv *.txt ./190524_alignment_summaries/

mv *.bam* ./190524_bam_files/

mv *.sam ./190524_sam_files/
############

echo "DONE!!"

##lastline
