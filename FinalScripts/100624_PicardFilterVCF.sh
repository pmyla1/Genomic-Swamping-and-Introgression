#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=2:00:00
#SBATCH --job-name=PicardFilterVCF
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


##load the Picard module for filtering the VCF
module load picard-uoneasy/3.0.0-Java-17

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/


###make environmental variables for the reference genome, input, and the output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
IN=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F1.biallelic.g.vcf.gz
OUT=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F2.best.practice.g.vcf.gz
################

##use Picard FilterVcf to hard-filter the vcf based on depth, quality, and allele balance
java -jar $EBROOTPICARD/picard.jar FilterVcf -R $REF -I $IN -O $OUT --MIN_AB 0.2 --MIN_DP 20 --MIN_GQ 25 --MIN_QD 2 --MAX_FS 60 
################

##unload the Picard module
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!!"

##lastline
