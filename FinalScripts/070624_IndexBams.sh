#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=12g
#SBATCH --time=01:00:00
#SBATCH --job-name=IndexBams
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load Samtools
module load samtools-uoneasy/1.18-GCC-12.3.0

############
meta=EKDL240001890-1A_222TKYLT4
##index duplicate marked bams with RG added
cd ~/220524_alignments/bam_files/duplicate_marked_bams_RG/

#first for FLE_2
samtools index -M -b -@ 4 ./FLEET_2_${meta}.marked_duplicates.bam 
#first for HAM_1
samtools index -M -b -@ 4 ./HAM_1_${meta}.marked_duplicates.bam
#first for NOT
samtools index -M -b -@ 4 ./NOT_${meta}.marked_duplicates.bam
#first for LWS
samtools index -M -b -@ 4 ./LWS_${meta}.marked_duplicates.bam 
#first for PAR_2
samtools index -M -b -@ 4 ./PAR_2_${meta}.marked_duplicates.bam
#first for Pen_1
samtools index -M -b -@ 4 ./Pen_1_${meta}.marked_duplicates.bam 
#first for SBAY_1
samtools index -M -b -@ 4 ./SBAY_1_${meta}.marked_duplicates.bam
#first for SPEY_2
samtools index -M -b -@ 4 ./SPEY_2_${meta}.marked_duplicates.bam 
#first for Iac
samtools index -M -b -@ 4 ./Iac.marked_duplicates.bam 
#first for Ime
samtools index -M -b -@ 4 ./Ime.marked_duplicates.bam

module unload samtools-uoneasy/1.18-GCC-12.3.0


echo "DONE!"

##lastline
