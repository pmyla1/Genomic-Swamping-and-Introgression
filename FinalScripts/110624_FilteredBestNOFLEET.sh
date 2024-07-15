#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=24:00:00
#SBATCH --job-name=Filtered.Best
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/

mkdir -p 110624_filtered.best/
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/300524_HaplotypeCaller_output/090624_combined_genotyped/090624_Ion.dan.NOFLEET2.genotyped.g.vcf.gz
OUT1=~/300524_HaplotypeCaller_output/090624_combined_genotyped/110624_filtered.best/110624_Ion.dan.F1.biallelic.g.vcf.gz
################
##GATK SelectVariants to select biallelic variants only
gatk SelectVariants \
   -R $REF \
   -V $VCF \
   -O $OUT1 \
   --select-type-to-exclude INDEL \
   --select-type-to-exclude MIXED \
   --restrict-alleles-to BIALLELIC 

##unload the GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##load the Picard module for filtering the VCF
module load picard-uoneasy/3.0.0-Java-17

cd ~/300524_HaplotypeCaller_output/090624_combined_genotyped/110624_filtered.best/

###make environmental variables for the reference genome, input, and the output
OUT2=~/300524_HaplotypeCaller_output/090624_combined_genotyped/110624_filtered.best/110624_Ion.dan.F2.best.g.vcf.gz
################

##use Picard FilterVcf to hard-filter the vcf based on depth, quality, and allele balance
java -jar $EBROOTPICARD/picard.jar FilterVcf -R $REF -I $OUT1 -O $OUT2 --MIN_AB 0.2 --MIN_DP 30 --MIN_GQ 30 --MIN_QD 2 --MAX_FS 60
################

##unload the Picard module
module unload picard-uoneasy/3.0.0-Java-17

##load vcftools
module load vcftools-uoneasy/0.1.16-GCC-12.3.0

#####
##use vcftools to output per site depth statistics on the hard-filtered (F2.best.practice) gVCF
vcftools --gzvcf $OUT2 --out 110624.depth.per.site --site-depth

module unload vcftools-uoneasy/0.1.16-GCC-12.3.0

echo "DONE!!"

##lastline
