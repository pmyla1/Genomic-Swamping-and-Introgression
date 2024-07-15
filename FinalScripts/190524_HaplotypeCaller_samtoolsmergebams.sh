#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=44g
#SBATCH --time=08:00:00
#SBATCH --job-name=HaplotypeCaller
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the GATK module for haplotype calling, specifying -ploidy 60 
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

########
##make a sequence dictionary for the C_excelsa_V5 reference genome
gatk CreateSequenceDictionary -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa

module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
########

#########
##load samtools to index fasta file
module load samtools-uoneasy/1.18-GCC-12.3.0

##index reference genome
samtools faidx ~/C_excelsa_V5_reference/C_excelsa_V5.fa

#######
##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/
##merge the duplicate_marked.bam files into one with samtools merge
samtools merge -@ 8 -o ./190524_extradanicas_ionopsidium_merged.sorted.bam ./190524_FLE_2_marked_duplicates.sorted.bam ./190524_HAM_1_marked_duplicates.sorted.bam ./190524_LWS_1_marked_duplicates.sorted.bam ./190524_NOT_5_marked_duplicates.sorted.bam ./190524_PAR_2_marked_duplicates.sorted.bam ./190524_PEN_1_marked_duplicates.sorted.bam ./190524_SBA_1_marked_duplicates.sorted.bam ./190524_SPE_2_marked_duplicates.sorted.bam ./190524_Iac_marked_duplicates.sorted.bam ./190524_Ime_marked_duplicates.sorted.bam 

module unload samtools-uoneasy/1.18-GCC-12.3.0
#########

#############
##re-load gatk module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

mkdir HaplotypeCaller_output/
#########
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_extradanicas_ionopsidium_merged.sorted.bam \
   -O ./HaplotypeCaller_output/190524_extra_danicas_Ionopsidium.g.vcf.gz \
   -ERC GVCF \
   --sample-ploidy 60 \
   --sites-only-vcf-output true
############

###unload module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
