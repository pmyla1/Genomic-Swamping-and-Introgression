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

###load the fastqc module for sequencing quality control
module load fastqc-uoneasy/0.12.1-Java-11
###############

################
##use fastqc specifying the output directory for the results
##Do this for the newly sequenced Cochlearia danica samples
##first on the FLE_2 sequencing data
fastqc -o ../170524_fastqc/ ./FLE_2/*.fq.gz
##next on HAM_1
fastqc -o ../170524_fastqc/ ./HAM_1/*.fq.gz
##next on LWS_1
fastqc -o ../170524_fastqc/ ./LWS_1/*.fq.gz
##next on NOT_5
fastqc -o ../170524_fastqc/ ./NOT_5/*.fq.gz
##next on PAR_2
fastqc -o ../170524_fastqc/ ./PAR_2/*.fq.gz
##next on PEN_1
fastqc -o ../170524_fastqc/ ./PEN_1/*.fq.gz
##next on SBA_1
fastqc -o ../170524_fastqc/ ./SBA_1/*.fq.gz
##next on SPE_2
fastqc -o ../170524_fastqc/ ./SPE_2/*.fq.gz

##now perform fastqc for the outgroup sequencing data
##Iac & Ime - Ionopsidium (sister taxa and outgroup for donwstream analyses - Dsuite)
fastqc -o ../170524_fastqc/ ./Iac/*.fastq.gz

fastqc -o ../170524_fastqc/ ./Ime/*.fastq.gz
##############

##unload fastqc module
module unload fastqc-uoneasy/0.12.1-Java-11
##########

#############
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a
##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc -f -p /gpfs01/home/pmyla1/170524_fastqc/.*fastqc.zip 

##unload multiqc module
module unload multiqc-uoneasy/1.14-foss-2023a

echo "DONE!!!"

##lastline