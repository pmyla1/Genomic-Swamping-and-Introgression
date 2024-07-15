#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=04:00:00
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

mkdir ~/220524_alignments/

##############
##try to make a for loop to iterate over the input files
#for infile in *_L1_1.trimmed.f*q.gz
#do
#      base=(${basename} ${infile} _L1_1.f*q.gz) 
#      bwa mem \
#           -t 8 ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
#           -R $RG \
#           ${infile} ${base}_L1_2.trimmed.f*q.gz \
#           > ~/220524_alignments/${base}_aln-pe.sam
##done

########################
##make a variable for storing the read group (RG) and read group ID (RGID)
#RGID=${sample}_${sample_id}_${library}_${flowcell}_${lane}_${barcode}
##read groups to add to the sam files for easier downstream analyses
#RG="@RG\tID:${RGID}\tSM:${sample}\tPL:${platform}\tLB:${library}\tPU:${flowcell}.${lane}.${barcode}"

#make a variable for the "metadata" 
metadata=EKDL240001890-1A_222TKYLT4

#store the output directory and reference genome in variables called OUT & REF, respectively
OUT=~/220524_alignments
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
###############
bwa mem \
     -t 8 $REF \
     -R $RG \
     ./Pen_1_${metadata}_L1_1.trimmed.fq.gz ./Pen_1_${metadata}_L1_2.trimmed.fq.gz \
     > $OUT/Pen_1_${metadata}_aln-pe.sam
##############
bwa mem \
     -t 8 $REF \
     -R $RG \
     ./SBAY_1_${metadata}_L1_1.trimmed.fq.gz ./SBAY_1_${metadata}_L1_2.trimmed.fq.gz \
     > $OUT/SBAY_1_${metadata}_aln-pe.sam
##############
bwa mem \
     -t 8 $REF \
     -R $RG \
     ./SPEY_2_${metadata}_L1_1.trimmed.fq.gz ./SPEY_2_${metadata}_L1_2.trimmed.fq.gz \
     > $OUT/SPEY_2_${metadata}_aln-pe.sam
#############
#bwa mem \
#     -t 8 $REF \
#     -R $RG \
#     ./Iac_1P.trimmed.fastq.gz ./Iac_2P.trimmed.fastq.gz \
#     > $OUT/Iac_aln-pe.sam
#############
bwa mem \
     -t 8 $REF \
     -R $RG \
     ./Ime_1P.trimmed.fastq.gz ./Ime_2P.trimmed.fastq.gz \
     > $OUT/Ime_aln-pe.sam
############

###unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

#############
##convert sam to bam, sort bams, index bams, get statistics
module unload samtools-uoneasy/1.18-GCC-12.3.0
###########
#on Pen_1
samtools view -@ 4 -h -b ./Pen_1_${metadata}_aln-pe.sam -o ./bam_files/Pen_1_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/Pen_1_${metadata}.sorted.bam ./bam_files/Pen_1_${metadata}.bam
samtools index ./bam_files/Pen_1_${metadata}.sorted.bam
samtools flagstat ./bam_files/Pen_1_${metadata}.sorted.bam > ./bam_files/Pen_1_${metadata}.flagstats
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
#on Iac
samtools view -@ 4 -h -b ./Iac_aln-pe.sam -o ./bam_files/Iac.bam
samtools sort -@ 4 -o ./bam_files/Iac.sorted.bam ./bam_files/Iac.bam
samtools index ./bam_files/Iac.sorted.bam
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
