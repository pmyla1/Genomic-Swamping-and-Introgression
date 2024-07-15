#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=MERGEVCFS
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load bcftools
module load bcftools-uoneasy/1.18-GCC-13.2.0

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/

###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF2OLD=~/reheadered.F4_133.ann.vcf.gz
VCF2NEW=~/110624_reheadered.F4_133.ann.vcf.gz
OUTMERGED=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_Merged.F4_133.ann.Ion.dan.filtered.F4.vcf.gz
VCF1OLD=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_Ion.dan.filtered.F4.g.vcf.gz
VCF1NEW=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/110624_reheadered.Ion.dan.filtered.F4.g.vcf.gz
#############
##use bcftools view to see the header, store in a new txt file called header{n}.txt so you can manually alter using nano
bcftools view -h $VCF1OLD > ./header1.txt

bcftools view -h $VCF2OLD > ./header2.txt
###########
##AT THIS POINT THE ORIGINAL SCRIPT ENDED HERE SO YOU CAN NAVIGATE TO THE DIRECTORY CONTAINING THE HEADER1.TXT AND HEADER2.TXT FILES AND MANUALLY ALTER BEFORE EXECUTING THE REST OF THE SCRIPT
#############
##manually altered the header1.txt and header2.txt PL fields to Number=.
bcftools reheader -h header1.txt $VCF1OLD > $VCF1NEW

bcftools reheader -h header2.txt $VCF2OLD > $VCF2NEW

#########
##index new VCF files
bcftools index $VCF1NEW

bcftools index $VCF2NEW
############
###########
##use bcftools merge, with 8 threads, assuming missing genotypes are 0/0, automatic indexing, and gzipped.
bcftools merge --t 8 -0 --write-index -Oz $VCF1NEW $VCF2NEW -o $OUTMERGED

##unload bcftools
module unload bcftools-uoneasy/1.18-GCC-13.2.0

echo "DONE!!"

##lastline
