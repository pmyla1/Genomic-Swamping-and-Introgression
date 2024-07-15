#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=8g
#SBATCH --time=01:00:00
#SBATCH --job-name=read_group_identification
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the samtools module to grep the read groups from the duplicate marked BAM files
module load samtools-uoneasy/1.18-GCC-12.3.0

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

#mkdir ~/2024.Cochlearia.Illumina.cohort/200524_read_groups/
########
##samtools view -H sample.bam | grep '^@RG' > ~/224.Cochlearia.Illumina.cohort/200524_read_groups/200524_all_bam_read.groups
#samtools view -H ./*duplicates.sorted.bam | grep '@RG' #> ~/2024.Cochlearia.Illumina.cohort/200524_read_groups/200524_all_bam_read_groups.txt

samtools view -H ./*duplicates.sorted.bam | grep '^@PG' > ~/2024.Cochlearia.Illumina.cohort/200524_read_groups/200524_all_bam_read_groups.txt


module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
