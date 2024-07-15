#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=32g
#SBATCH --time=06:00:00
#SBATCH --job-name=BWA_align_reads
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##cd to C_excelsa_V5_reference
#cd /gpfs01/home/pmyla1/C_excelsa_V5_reference/

#######
##load bwa module
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

########
##use BWA to index the C_excelsa_V5.fa reference
#bwa index ./C_excelsa_V5.fa

##use bwa mem to align the short reads from the FLE_2 population to the C_excelsa_V5.fa 
#bwa mem ./C_excelsa_V5.fa ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ../2024.Cochlearia.Illumina.cohort/170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz > ../160524_new_danica_alignments/160524_aln-se_FLE_2.sam


#module load bwa-uoneasy/0.7.17-GCCcore-12.3.0

##cd to 170524_cutadapt/
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/170524_cutadapt/180524_merged_reads/

##make a new directory to store the alignments
#mkdir 180524_alignments/

##make a test directory for the first alignment
#mkdir 180524_test_alignment/
##use bwa mem on one of the merged fastq.gz files to align to the reference
##use the -p flag for "interleaved" reads (i.e. R1 and R2 paired end reads merged into a single file)
bwa mem -t 16 /gpfs01/home/pmyla1/C_excelsa_V5_reference/C_excelsa_V5.fa ./180524_merged_FLEET_2_EKDL240001890-1A_222TKYLT4.fq.gz > ./180524_test_alignment/180524_FLEET_2_EKDL240001890-1A_222TKYLT4_paired.sam


#module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0
##unload bwa module
module unload bwa-uoneasy/0.7.17-GCCcore-12.3.0

echo "DONE!!"

##lastline
