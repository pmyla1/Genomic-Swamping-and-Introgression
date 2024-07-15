#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=Samtools_bam_stats
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the picard module to mark the duplicates in the bam alignment files
module load samtools-uoneasy/1.18-GCC-12.3.0

##change directory to the bam alignment files
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

mkdir 190524_samtools_stats
#########
##use samtools-stats to provide metrics on the duplicate-marked bam alignments
##first on FLE_2
samtools stats ./190524_FLE_2_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_FLE_2.bam.stats
##next on HAM_1
samtools stats ./190524_HAM_1_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_HAM_1.bam.stats
##then on LWS_1 
samtools stats ./190524_LWS_1_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_LWS_1.bam.stats
##now on NOT_5
samtools stats ./190524_NOT_5_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_NOT_5.bam.stats
##next on PAR_2
samtools stats ./190524_PAR_2_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_PAR_2.bam.stats
##then on PEN_1
samtools stats ./190524_PEN_1_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_PEN_1.bam.stats
##now on SBA_1
samtools stats ./190524_SBA_1_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_SBA_1.bam.stats
##finally on SPE_2
samtools stats ./190524_SPE_2_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_SPE_2.bam.stats

############
##now do the same for the Ionopsidium samples
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/Ionopsidium_cutadapt/180524_Ionopsidium_merged_reads/180524_Ionopsidium_alignments/190524_bam_files/

mkdir 190524_samtools_stats

##first on Iac
samtools stats ./190524_Iac_marked_duplicates.sorted.bam > ./190524_samtools_stats/190524_Iac.bam.stats
##finally on Ime
samtools stats ./190524_Ime_marked_duplicates.sorted.bam  > ./190524_samtools_stats/190524_Ime.bam.stats

#########

###unload module
module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
