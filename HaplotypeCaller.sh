#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=HaplotypeCaller
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the GATK module for haplotype calling, specifying -ploidy 6 
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

########
##make a sequence dictionary for the C_excelsa_V5 reference genome
gatk CreateSequenceDictionary -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa

module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

#########
##load samtools to index fasta file
module load samtools-uoneasy/1.18-GCC-12.3.0

samtools faidx ~/C_excelsa_V5_reference/C_excelsa_V5.fa

module unload samtools-uoneasy/1.18-GCC-12.3.0
#########

#############
##re-load gatk module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##change directory to the bam alignment files
cd ~/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/180524_alignments/190524_bam_files/

#mkdir HaplotypeCaller_output/
#########
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_FLE_2_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.FLE_2.g.vcf.gz
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_HAM_1_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.HAM_1.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6 
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_LWS_1_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.LWS_1.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_NOT_5_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.NOT_5.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_PAR_2_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.PAR_2.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_LWS_1_marked_duplicates.sorted.bam 
   -O ./HaplotypeCaller_output/190524.LWS_1.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_PEN_1_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.PEN_1.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
###########
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_SBA_1_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.SBA_1.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_SPE_2_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.SPE_2.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6 
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_Iac_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.Iac.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6 
############
##use GATK HaplotypeCaller to make a gVCF of the additional C. danica and Ionopsidium files
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ~/C_excelsa_V5_reference/C_excelsa_V5.fa \
   -I ./190524_Ime_marked_duplicates.sorted.bam \
   -O ./HaplotypeCaller_output/190524.Ime.g.vcf.gz -ERC GVCF \
   --sample-ploidy 6 

############

###unload module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline
