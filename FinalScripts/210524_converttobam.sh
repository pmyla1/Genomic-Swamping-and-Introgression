#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=24g
#SBATCH --time=04:00:00
#SBATCH --job-name=converttobam
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile
#######

##change directory to the 180524_alignments
cd ~/200524_new_alignments/
##load samtools view to convert the sam files to bam files
module load samtools-uoneasy/1.18-GCC-12.3.0

##########
##make directories to store the output 
mkdir 210524_sam_files/

mkdir 210524_bam_files/

mkdir 210524_alignment_summaries/
###########

##############
##convert FLE_2 sam to bam 
samtools view -@ 8 -b ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.sam > ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.bam
##use samtools sort to produce a sorted FLE_2 bam file
samtools sort -@ 8 -o ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_PE.sorted.bam ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.bam
##finally index the sorted FLE_2 bam file with samtools index
samtools index ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_PE.sorted.bam
##now get a summary of the FLE_2 alignment with samtools flagstat
samtools flagstat ./200524_FLEET_2_EKDL240001890-1A_222TKYLT4_PE.sorted.bam > ./200524_FLEET_2_aln_summary.txt


##convert HAM_1 sam to bam 
samtools view -@ 8 -b ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_paired.sam > ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_paired.bam
##use samtools sort to produce a sorted HAM_1 bam file
samtools sort -@ 8 -o ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_PE.sorted.bam ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_paired.bam
##finally index the sorted HAM_1 bam file with samtools index
samtools index ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_PE.sorted.bam
##now get a summary of the HAM_1 alignment with samtools flagstat
samtools flagstat ./200524_HAM_1_EKDL240001890-1A_222TKYLT4_PE.sorted.bam > ./200524_HAM_1_aln_summary.txt

##convert LWS_1 sam to bam 
samtools view -@ 8 -b ./200524_LWS_EKDL240001890-1A_222TKYLT4_paired.sam > ./200524_LWS_EKDL240001890-1A_222TKYLT4_paired.bam
##use samtools sort to produce a sorted LWS_1 bam file
samtools sort -@ 8 -o ./200524_LWS_EKDL240001890-1A_222TKYLT4_PE.sorted.bam ./200524_LWS_EKDL240001890-1A_222TKYLT4_paired.bam
##finally index the sorted LWS_1 bam file with samtools index
samtools index ./200524_LWS_EKDL240001890-1A_222TKYLT4_PE.sorted.bam
##now get a summary of the HAM_1 alignment with samtools flagstat
samtools flagstat ./200524_LWS_EKDL240001890-1A_222TKYLT4_PE.sorted.bam > ./200524_LWS_1_aln_summary.txt

##convert NOT_5 sam to bam
#samtools view -@ 8 -b ./180524_NOT_EKDL240001890-1A_222TKYLT4_paired.sam > ./180524_NOT_5_paired.bam
##use samtools sort to produce a sorted NOT_5 bam file
#samtools sort -@ 8 -o ./180524_NOT_5_paired.sorted.bam ./180524_NOT_5_paired.bam                             
##finally index the sorted NOT_5 bam file with samtools index
#samtools index ./180524_HAM_1_paired.sorted.bam
##now get a summary of the NOT_5 alignment with samtools flagstat
#samtools flagstat ./180524_NOT_5_paired.sorted.bam > ./180524_NOT_5_alignment_summary.txt

##convert PAR_2 sam to bam 
#samtools view -@ 8 -b ./180524_PAR_2_EKDL240001890-1A_222TKYLT4_paired.sam > ./180524_PAR_2_paired.bam
##use samtools sort to produce a sorted PAR_2 bam file
#samtools sort -@ 8 -o ./180524_PAR_2_paired.sorted.bam ./180524_PAR_2_paired.bam                             
##finally index the sorted PAR_2 bam file with samtools index
#samtools index ./180524_PAR_2_paired.sorted.bam
##now get a summary of the PAR_2 alignment with samtools flagstat
#samtools flagstat ./180524_PAR_2_paired.sorted.bam > ./180524_PAR_2_alignment_summary.txt

##convert PEN_1 sam to bam 
#samtools view -@ 8 -b ./180524_Pen_1_EKDL240001890-1A_222TKYLT4_paired.sam > ./180524_PEN_1_paired.bam
##use samtools sort to produce a sorted PEN_1 bam file
#samtools sort -@ 8 -o ./180524_PEN_1_paired.sorted.bam ./180524_PEN_1_paired.bam                             
##finally index the sorted PEN_1 bam file with samtools index
#samtools index ./180524_PEN_1_paired.sorted.bam
##now get a summary of the PEN_1 alignment with samtools flagstat
#samtools flagstat ./180524_PEN_1_paired.sorted.bam > ./180524_PEN_1_alignment_summary.txt

##convert SBA_1 sam to bam 
#samtools view -@ 8 -b ./180524_SBAY_1_EKDL240001890-1A_222TKYLT4_paired.sam > ./180524_SBA_1_paired.bam
##use samtools sort to produce a sorted SBA_1 bam file
#samtools sort -@ 8 -o ./180524_SBA_1_paired.sorted.bam ./180524_SBA_1_paired.bam                             
##finally index the sorted SBA_1 bam file with samtools index
#samtools index ./180524_SBA_1_paired.sorted.bam
##now get a summary of the SBA_1 alignment with samtools flagstat
#samtools flagstat ./180524_SBA_1_paired.sorted.bam > ./180524_SBA_1_alignment_summary.txt

##convert SPE_2 sam to bam 
#samtools view -@ 8 -b ./180524_SPEY_2_EKDL240001890-1A_222TKYLT4_paired.sam > ./180524_SPE_2_paired.bam
##use samtools sort to produce a sorted SPE_2 bam file
#samtools sort -@ 8 -o ./180524_SPE_2_paired.sorted.bam ./180524_SPE_2_paired.bam
##finally index the sorted SPE_2 bam file with samtools index
#samtools index ./180524_SPE_2_paired.sorted.bam
##now get a summary of the SPE_2 alignment with samtools flagstat
#samtools flagstat ./180524_SPE_2_paired.sorted.bam > ./180524_SPE_2_alignment_summary.txt

##convert Ime sam to bam
samtools view -@ 8 -b ./200524_Iac_paired.sam > ./200524_Iac_paired.bam
##sort the bam file
samtools sort -@ 8 -o ./200524_Iac_PE.sorted.bam ./200524_Iac_paired.bam
##index
samtools index ./200524_Iac_PE.sorted.bam
##flatstat
samtools flagstat ./200524_Iac_PE.sorted.bam > ./200524_Iac_aln_summary.txt

module unload samtools-uoneasy/1.18-GCC-12.3.0
###########

#############
##move the output files into their appropriate directories
#mv *.txt ./190524_alignment_summaries/

#mv *.sam ./190524_sam_files/

#mv *.bam* ./190524_bam_files/


echo "DONE!!"

##lastline
