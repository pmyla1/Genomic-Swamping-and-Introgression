#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
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

mkdir ~/220524_alignments/

########################
#make a variable for the "metadata" 
metadata=EKDL240001890-1A_222TKYLT4

##store the output directory and reference genome in variables called OUT & REF, respectively
OUT=~/220524_alignments
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa

####################
##bwa mem
bwa mem \
     -t 16 $REF \
     ./FLEET_2_*_L1_1.trimmed.fq.gz ./FLEET_2_*_L1_2.trimmed.fq.gz \
     > $OUT/FLEET_2_EKDL240001890-1A_222TKYLT4_aln-pe.sam
#################
bwa mem \
     -t 16 $REF \
     ./HAM_1_*_L1_1.trimmed.fq.gz ./HAM_1_*_L1_2.trimmed.fq.gz \
     > $OUT/HAM_1_EKDL240001890-1A_222TKYLT4_aln-pe.sam
###############
bwa mem \
     -t 16 $REF \
     ./LWS_*_L1_1.trimmed.fq.gz ./LWS_*_L1_2.trimmed.fq.gz \
     > $OUT/LWS_EKDL240001890-1A_222TKYLT4_aln-pe.sam
###############
bwa mem \
     -t 16 $REF \
     ./NOT_*_L1_1.trimmed.fq.gz ./NOT_*_L1_2.trimmed.fq.gz \
     > $OUT/NOT_EKDL240001890-1A_222TKYLT4_aln-pe.sam
###############
bwa mem \
     -t 16 $REF \
     ./PAR_2_${metadata}_L1_1.trimmed.fq.gz ./PAR_2_${metadata}_L1_2.trimmed.fq.gz \
     > $OUT/PAR_2_${metadata}_aln-pe.sam
###############
bwa mem \
     -t 16 $REF \
     ./Pen_1_*_L1_1.trimmed.fq.gz ./Pen_1_*_L1_2.trimmed.fq.gz \
     > $OUT/Pen_1_EKDL240001890-1A_222TKYLT4_aln-pe.sam
##############
bwa mem \
     -t 16 $REF \
     ./SBAY_1_*_L1_1.trimmed.fq.gz ./SBAY_1_*_L1_2.trimmed.fq.gz \
     > $OUT/SBAY_1_EKDL240001890-1A_222TKYLT4_aln-pe.sam
##############
bwa mem \
     -t 16 $REF \
     ./SPEY_2_${metadata}_L1_1.trimmed.fq.gz ./SPEY_2_${metadata}_L1_2.trimmed.fq.gz \
     > $OUT/SPEY_2_${metadata}_aln-pe.sam
#############
bwa mem \
     -t 16 $REF \
     -R $RG \
     ./Iac_1P.trimmed.fastq.gz ./Iac_2P.trimmed.fastq.gz \
     > $OUT/Iac_aln-pe.sam
#############
bwa mem \
     -t 16 $REF \
     -R $RG \
     ./Ime_1P.trimmed.fastq.gz ./Ime_2P.trimmed.fastq.gz \
     > $OUT/Ime_aln-pe.sam
############

###unload module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

#################
##CONVERT SAM FILES TO BAM FILES, SORT BAMS, INDEX, THEN FLAGSTAT WITH SAMTOOLS
module load samtools-uoneasy/1.18-GCC-12.3.0

##change directory to the 220524_alignments folder
cd ~/220524_aligments/

mkdir bam_files/

mkdir bam_files/sorted_bams/

mkdir bam_files/sorted_bams/samtools_flagstats/

##########
##samtools view with 4 threads, output headers(-h), in bam format (-b)
##samtools sort with 4 threads, using .bam as input
##samtool index
##samtools flagstat on the sorted bam
################
#on FLE_2
samtools view -@ 4 -h -b ./FLEET_2_${metadata}_aln-pe.sam -o ./bam_files/FLEET_2_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/FLEET_2_${metadata}.sorted.bam ./bam_files/FLEET_2_${metadata}.bam
samtools index ./bam_files/FLEET_2_${metadata}.sorted.bam
samtools flagstat ./bam_files/FLEET_2_${metadata}.sorted.bam > ./bam_files/FLEET_2_${metadata}.flagstats
########
##for HAM_1
samtools view -@ 4 -h -b ./HAM_1_${metadata}_aln-pe.sam -o ./HAM_1_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/HAM_1_${metadata}.sorted.bam ./HAM_1_${metadata}.bam
samtools index ./bam_files/HAM_1_${metadata}.sorted.bam
samtools flagstat ./bam_files/HAM_1_${metadata}.sorted.bam > ./bam_files/HAM_1_${metadata}.flagstats
###########
#on LWS_1
samtools view -@ 4 -h -b ./LWS_${metadata}_aln-pe.sam -o ./bam_files/LWS_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/LWS_${metadata}.sorted.bam ./bam_files/LWS_${metadata}.bam
samtools index ./bam_files/LWS_${metadata}.sorted.bam
samtools flagstat ./bam_files/LWS_${metadata}.sorted.bam > ./bam_files/LWS_${metadata}.flagstats
############
#on NOT_5
samtools view -@ 4 -h -b ./NOT_${metadata}_aln-pe.sam -o ./bam_files/NOT_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/NOT_${metadata}.sorted.bam ./bam_files/NOT_${metadata}.bam
samtools index ./bam_files/NOT_${metadata}.sorted.bam
samtools flagstat ./bam_files/NOT_${metadata}.sorted.bam > ./bam_files/NOT_${metadata}.flagstats
###########
#on PAR_2
samtools view -@ 4 -h -b ./PAR_2_${metadata}_aln-pe.sam -o ./bam_files/PAR_2_${metadata}.bam 
samtools sort -@ 4 -o ./bam_files/PAR_2_${metadata}.sorted.bam ./bam_files/PAR_2_${metadata}.bam
samtools index ./bam_files/PAR_2_${metadata}.sorted.bam
samtools flagstat ./bam_files/PAR_2_${metadata}.sorted.bam > ./bam_files/PAR_2_${metadata}.flagstats
###########
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
samtools sort -@ 4 -o ./bam_files/sorted_bams/Iac.sorted.bam ./bam_files/Iac.bam
samtools index ./bam_files/sorted_bams/Iac.sorted.bam
samtools flagstat ./bam_files/sorted_bams/Iac.sorted.bam > ./bam_files/sorted_bams/samtools_flagstats/Iac.flagstats
##########
#on Ime
samtools view -@ 4 -h -b ./Ime_aln-pe.sam -o ./bam_files/Ime.bam
samtools sort -@ 4 -o ./bam_files/sorted_bams/Ime.sorted.bam ./bam_files/Ime.bam
samtools index ./bam_files/sorted_bams/Ime.sorted.bam
samtools flagstat ./bam_files/sorted_bams/Ime.sorted.bam > ./bam_files/sorted_bams/samtools_flagstats/Ime.flagstats
###########
module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
