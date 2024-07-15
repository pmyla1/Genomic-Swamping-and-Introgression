#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=FASTQC
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##cd to the 2024.Cochlearia.Illumina.cohort folder where the fastq.gz files are 
cd /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/

mkdir Ionopsidium_trimmed_fastqc/

mkdir Ionopsidium_cutadapt/
##############
##cut adapters from the end of the fastq.gz files with cutadapt
##the Illumina universal adapter sequence is CTGTCTCTTATACACATCT
##the Illumina universal adapter sequence is CTGTCTCTTATACACATCT
module load cutadapt-uon/gcc12.3.0/4.6

##first cut the Illumina universal adapter sequence from Iac1P
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Iac_1P.fastq.gz ./Iac/Iac_1P.fastq.gz
##then cut the Illumina universal adapter sequence from Iac2P
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Iac_2P.fastq.gz ./Iac/Iac_2P.fastq.gz

##cut the Illumina universal adapter seqeucne from Iac1U
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Iac_1U.fastq.gz ./Iac/Iac_1U.fastq.gz
##then cut Illumina universal adapter sequence from Iac2U
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Iac_2U.fastq.gz ./Iac/Iac_2U.fastq.gz

##cut the Illumina universal adapter sequence from Ime1P
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Ime_1P.fastq.gz ./Ime/Ime_1P.fastq.gz
##then cut the Illumina universal adapter sequence from Ime2p
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Ime_2P.fastq.gz ./Ime/Ime_2P.fastq.gz

##cut the Illumina universal adapter sequence from Ime1U
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Ime_1U.fastq.gz ./Ime/Ime_1U.fastq.gz
##then cut the Illumina universal adapter sequence from Ime2U
cutadapt -a CTGTCTCTTATACACATCT -o Ionopsidium_cutadapt/170524_Ime_2U.fastq.gz ./Ime/Ime_2U.fastq.gz

##unload cutadapt module
module unload cutadapt-uon/gcc12.3.0/4.6

#############
##SEQUENCING QUALITY CONTROL
###load the fastqc module for sequencing quality control
module load fastqc-uoneasy/0.12.1-Java-11
###############

################
##use fastqc specifying the output directory for the results
##Do this for the newly sequenced Cochlearia danica samples
##specify the output to go into the 170524_fastqc_trimmed_output folder, and to execute fastqc on the 170524_cutadapt/ directory
fastqc -o ./Ionopsidium_trimmed_fastqc/ ./Ionopsidium_cutadapt/*fastq.gz
#############

##unload fastqc module
module unload fastqc-uoneasy/0.12.1-Java-11
############

###############
##MULTIPLE QUALITY CONTROL REPORTS ON THE FASTQC.ZIP FILES PRODUCED BY FASTQC
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a
##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc -f -p /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/Ionopsidium_trimmed_fastqc/*fastqc.zip 

##unload multiqc module
module unload multiqc-uoneasy/1.14-foss-2023a
###############

echo "DONE!!!"

##lastline
