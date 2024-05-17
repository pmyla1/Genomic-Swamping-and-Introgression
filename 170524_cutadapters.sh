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

##make a new directory for the cutadapt output 
mkdir 170524_cutadapt/

##make a new directory for the trimmed fastqc output
mkdir 170524_fastqc_trimmed_output/
##############
##cut adapters from the end of the fastq.gz files with cutadapt
##the read 1 nextera transposase sequence is TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
## the read 2 nextera transposase adapter sequence is GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
module load cutadapt-uon/gcc12.3.0/4.6

##first cut the nextera transposase read 1 sequence from FLE_2
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from FLE_2
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz

##cut the nextera transposase read 1 sequence from HAM_1
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from HAM_1
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from LWS_1
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./LWS_1/LWS_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from LWS_1
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./LWS_1/LWS_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from NOT_5
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./NOT_5/NOT_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from NOT_5
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./NOT_5/NOT_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from PAR_2
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./PAR_2/PAR_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from PAR_2
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./PAR_2/PAR_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from PEN_1
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 1170524_cutadapt/170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./PEN_1/Pen_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from PEN_1
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./PEN_1Pen_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from SBA_1
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./SBA_1/SBAY_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from SBA_2
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./SBA_1/SBAY_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

##cut the nextera transposase read 1 sequence from SPE_2
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./SPE_2/SPEY_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz
##then cut the nextera transposase read 2 sequence from SPE_2
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./SPE_2/SPEY_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz

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
fastqc -o ../170524_fastqc_trimmed_output/ ./170524_cutadapt/*.fq.gz
#############

##unload fastqc module
module unload fastqc-uoneasy/0.12.1-Java-11
############

###############
##MULTIPLE QUALITY CONTROL REPORTS ON THE FASTQC.ZIP FILES PRODUCED BY FASTQC
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a
##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc -f -p /gpfs01/home/pmyla1/2024.Cochlearia.Illumina.cohort/170524_fastqc_trimmed_output/.*fastqc.zip 

##unload multiqc module
module unload multiqc-uoneasy/1.14-foss-2023a
###############

echo "DONE!!!"

##lastline