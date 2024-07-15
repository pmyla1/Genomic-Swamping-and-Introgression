#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=08:00:00
#SBATCH --job-name=MAPREADSSAMTOOLSSORT
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the trimmomatic module 
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0


##change directory to 220524_trimmed_reads
cd ~/2024.Cochlearia.Illumina.cohort/220524_trimmed_reads/
############
#make a variable for the "metadata" 
metadata=EKDL240001890-1A_222TKYLT4

#store the output directory and reference genome in variables called OUT & REF, respectively
OUT=~/220524_alignments
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa

####################

#############
##Iac first
#bwa mem \
#     -t 16 $REF \
#     ./Iac_1P.trimmed.fastq.gz ./Iac_2P.trimmed.fastq.gz \
#     > $OUT/Iac_aln-pe.sam
###########
##Ime last
bwa mem \
     -t 16 $REF \
     ./Ime_1P.trimmed.fastq.gz ./Ime_2P.trimmed.fastq.gz \
     > $OUT/Ime_aln-pe.sam
############

###unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0


#############
##samtools to convert to bam, sort, index, and get alignment statistics
module load samtools-uoneasy/1.18-GCC-12.3.0

cd ~/220524_alignments/
#on Iac
#samtools view -@ 4 -h -b ./Iac_aln-pe.sam -o ./bam_files/Iac.bam
#samtools sort -@ 4 -o ./bam_files/Iac.sorted.bam ./bam_files/Iac.bam
#samtools index ./bam_files/Iac.sorted.bam
samtools flagstat ./bam_files/Iac.sorted.bam > ./bam_files/Iac.flagstats
##########
#on Ime
samtools view -@ 4 -h -b ./Ime_aln-pe.sam -o ./bam_files/Ime.bam
samtools sort -@ 4 -o ./bam_files/Ime.sorted.bam ./bam_files/Ime.bam
samtools index ./bam_files/Ime.sorted.bam
samtools flagstat ./bam_files/Ime.sorted.bam > ./bam_files/Ime.flagstats
###########
module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
