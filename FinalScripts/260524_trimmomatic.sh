#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=40g
#SBATCH --time=06:00:00
#SBATCH --job-name=Trimmomatic
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

############
##load the trimmomatic module 
module load trimmomatic-uoneasy/0.39-Java-17

##make a new directory to store the trimmed reads 
mkdir ~/2024.Cochlearia.Illumina.cohort/220524_trimmed_reads/

##change directory to 220524_trimmed_reads
cd ~/2024.Cochlearia.Illumina.cohort/220524_trimmed_reads/

###################
##use trimmomatic on paired-end (PE) reads, with phred33 scoring, and nextera transposase PE-adapters
##using a sliding window of 4bp and removing reads with phred scores < 20 
########################
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../FLE_2/*_L1_1.fq.gz ../FLE_2/*_L1_2.fq.gz \
        ./FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
##############
##cycle through the different accessions
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../HAM_1/*_L1_1.fq.gz ../HAM_1/*_L1_2.fq.gz \
        ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
###############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../LWS_1/*_L1_1.fq.gz ../LWS_1/*_L1_2.fq.gz \
        ./LWS_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./LWS_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./LWS_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./LWS_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
###############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../NOT_5/*_L1_1.fq.gz ../NOT_5/*_L1_2.fq.gz \
        ./NOT_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./NOT_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./NOT_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./NOT_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
###############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../PAR_2/*_L1_1.fq.gz ../PAR_2/*_L1_2.fq.gz \
        ./PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
###############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../PEN_1/*_L1_1.fq.gz ../PEN_1/*_L1_2.fq.gz \
        ./Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
##############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../SBA_1/*_L1_1.fq.gz ../SBA_1/*_L1_2.fq.gz \
        ./SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
##############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../SPE_2/*_L1_1.fq.gz ../SPE_2/*_L1_2.fq.gz \
        ./SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
#############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../Iac/Iac_1P.fastq.gz ../Iac/Iac_2P.fastq.gz \
        ./Iac_1P.trimmed.fastq.gz ./Iac_1P.un.orhpan.fastq.gz \
        ./Iac_2P.trimmed.fastq.gz ./Iac_2P.un.orhpan.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
#############
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../Ime/Ime_1P.fastq.gz ../Ime/Ime_2P.fastq.gz \
        ./Ime_1P.trimmed.fastq.gz ./Ime_1P.orhpan.fastq.gz \
        ./Ime_2P.trimmed.fastq.gz ./Ime_2P.orphan.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
#############


###unload module
module unload trimmomatic-uoneasy/0.39-Java-17

echo "DONE!!"

##lastline
