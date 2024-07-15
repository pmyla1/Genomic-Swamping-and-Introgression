#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --job-name=DepthOfCoverage
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

####load samtools to index duplicate marked bam files
module load samtools-uoneasy/1.18-GCC-12.3.0

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

#################
##index the duplicate-marked files with samtools index
samtools index ./190524_FLE_2_marked_duplicates.sorted.bam
samtools index ./190524_HAM_1_marked_duplicates.sorted.bam
samtools index ./190524_LWS_1_marked_duplicates.sorted.bam
samtools index ./190524_NOT_5_marked_duplicates.sorted.bam
samtools index ./190524_PAR_2_marked_duplicates.sorted.bam
samtools index ./190524_PEN_1_marked_duplicates.sorted.bam
samtools index ./190524_SBA_1_marked_duplicates.sorted.bam
samtools index ./190524_SPE_2_marked_duplicates.sorted.bam
samtools index ./190524_Iac_marked_duplicates.sorted.bam
samtools index ./190524_Ime_marked_duplicates.sorted.bam

module unload samtools-uoneasy/1.18-GCC-12.3.0
############
##load the GATK module to grep the read groups from the duplicate marked BAM files
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

#############
##calculate Depth of Coverage for the Input bam files (duplicate marked)
 gatk \
   DepthOfCoverage \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -O 200524_DepthOfCoverage \
   -I ./200524_input_bams.list \
   -L ~/2024.Cochlearia.Illumina.cohort/200524_Intervals.list
#############

module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
