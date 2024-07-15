#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --job-name=MergeVCFs
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load Picard
module load picard-uoneasy/3.0.0-Java-17

cd ~/280524_final_vcfs/

mkdir ~/280524_Ionops_Cochlearia_output/

OUTDIR=~/280524_Ionops_Cochlearia_output
############
##merge the 4dg_ionops_only_filtered.F4.vcf.gz (provided by Yant, 28/05/2024) with the reheadered.F4_133.ann.vcf.gz (also provided by Yant, 09/05/2024)
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
          I=./4dg_ionops_only_filtered.F4.vcf.gz \
          I=~/reheadered.F4_133.ann.vcf.gz \
          O=$OUTDIR/280524_Ionops_Cochlearia_combined.vcf.gz
###########

##unload Picard
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline
