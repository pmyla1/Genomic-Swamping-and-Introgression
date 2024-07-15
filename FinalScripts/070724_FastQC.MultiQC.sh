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

mkdir -p ~/070724.FastQC.MultiQC
##cd to where the fastq.gz files for the outgroup (Ionopsidium) samples are 
cd ~/ionops/01.merged/

###load the fastqc module for sequencing quality control
module load fastqc-uoneasy/0.12.1-Java-11
###############

################
##use fastqc specifying the output directory for the results
##firstly perform FastQC for the Iab_1 and Iab_2 data
fastqc -o ~/070724.FastQC.MultiQC/ ./Iab_1/*.fastq.gz

fastqc -o ~/070724.FastQC.MultiQC/ ./Iab_2/*.fastq.gz

fastqc -o ~/070724.FastQC.MultiQC/ ./Iac/*.fastq.gz

fastqc -o ~/070724.FastQC.MultiQC/ ./Ime/*.fastq.gz
#############
##change directory to the additional Danica samples
cd ~/2024.Cochlearia.Illumina.cohort/

##Do this for the newly sequenced Cochlearia danica samples
##first on the FLE_2 sequencing data
fastqc -o ~/070724.FastQC.MultiQC/ ./FLE_2/*.fq.gz
##next on HAM_1
fastqc -o ~/070724.FastQC.MultiQC/ ./HAM_1/*.fq.gz
##next on LWS_1
fastqc -o ~/070724.FastQC.MultiQC/ ./LWS_1/*.fq.gz
##next on NOT_5
fastqc -o ~/070724.FastQC.MultiQC/ ./NOT_5/*.fq.gz
##next on PAR_2
fastqc -o ~/070724.FastQC.MultiQC/ ./PAR_2/*.fq.gz
##next on PEN_1
fastqc -o ~/070724.FastQC.MultiQC/ ./PEN_1/*.fq.gz
##next on SBA_1
fastqc -o ~/070724.FastQC.MultiQC/ ./SBA_1/*.fq.gz
##next on SPE_2
fastqc -o ~/070724.FastQC.MultiQC/ ./SPE_2/*.fq.gz
##################
##unload fastqc module
module unload fastqc-uoneasy/0.12.1-Java-11
##########

#################
##multiqc analysis of C. danica sequencing data
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a
##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc -f -p ~/070724.FastQC.MultiQC/*fastqc.zip  

##unload multiqc module
module unload multiqc-uoneasy/1.14-foss-2023a


echo "DONE!!!"

##lastline
