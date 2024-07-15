#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=2:00:00
#SBATCH --job-name=MeanDepthStats
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/
##load python module
module load python-uoneasy/3.11.5-GCCcore-13.2.0

###make environmental variables for the reference genome, input, and the output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
IN=~/300524_HaplotypeCaller_output/090624_combined_genotyped/Depth.per.site.ldepth
OUT=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_meandepthstats.txt

##########
import statistics 
input=open($IN, "r")
output=open($OUT,"w") 
depths=[int(line.split('\t')[2]) for line in input]
mean_depth=statistics.mean(depths) ##calculate the mean depth
cut_off=int(round(mean_depth)*1.6) #to calculate the upper limit for depth and round up to match the DP=int field 
##close input and output files
input.close()
output.close()
##unload python module
module unload python-uoneasy/3.11.5-GCCcore-13.2.0

echo "DONE!!"

##lastline
