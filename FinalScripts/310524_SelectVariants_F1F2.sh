#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40g
#SBATCH --time=08:00:00
#SBATCH --job-name=SelectVariantsFilter
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile


##load the GATK module for combining the GVCFs
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

cd 310524_Genotyped_VCF/

mkdir -p 310524_filtered_VCFs
###make environmental variables for the reference genome, the input VCF and the F1 output VCF
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/310524_Genotyped_VCF/310524_Genotyped_Ionopsidium.vcf.gz
F1_OUT=~/310524_Genotyped_VCF/310524_filtered_F1_biallelic_Ionopsidium.vcf.gz
F2_OUT=~/310524_Genotyped_VCF/310524_filtered_F2_bestpractice_Ionopsidium.vcf.gz
################
##GATK SelectVariants - exclude insertion/deletion mutations, exclude variants with mixed SNPs/indels at the same site
##restrict alleles to biallelic SNPs.
gatk SelectVariants \
 -R $REF \
 -V $VCF \
 -O $F1_OUT \
 --select-type-to-exclude INDEL \
 --select-type-to-exclude MIXED \
 --restrict-alleles-to BIALLELIC
##############   

##############
##GATK VariantFiltration
##filter with a quality score normalised by unfiltered allele depth (QUAL/DP) < 2.0 
##FS is Fisher's exact test to detect forward/reverse strand allele bias = artifacts due to sequencing errors
##MQ is the root mean squared mapping quality
##MQRS detects read mapping quality allele bias
##RPRS detect relative read position allele bias
##HS filters when multiple segregating haplotypes
gatk VariantFiltration \
    -R $REF \
    -V $F1_OUT \
    -O $F2_OUT \
    --filter-name "QD" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS" \
    --filter-expression "FS > 60" \
    --filter-name "MQ" \
    --filter-expression "MQ < 40" \
    --filter-name "MQRS" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "RPRS" \
    --filter-expression "ReadPosRankSum < -8" \
    --filter-name "HS" \
    --filter-expression "HaplotypeScore < 13"
#############

##unload the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

#############
##load vcftools
module load vcftools-uoneasy/0.1.16-GCC-12.3.0

###########
##use vcftools to output depth-per-site statistics on the F2 filtered vcf
vcftools --gzvcf $F2_OUT --out ~/310524_Genotyped_VCF/310524_F2_DepthPerSite --site-depth

module unload vcftools-uoneasy/0.1.16-GCC-12.3.0
##########

echo "DONE!!"

##lastline
