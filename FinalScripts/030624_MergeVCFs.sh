#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100g
#SBATCH --time=01:00:00
#SBATCH --job-name=MergeVCFs
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load bcftools
#module load bcftools-uoneasy/1.18-GCC-13.2.0

############
##make environmental variables for the input VCFs to be sorted, index, and merged
#INVCF1=~/reheadered.F4_133.ann.vcf.gz
#OUTVCF1=~/030624_mergedVCF/sorted_reheadered.F4_133.ann.vcf.gz

#INVCF2=~/ionops/10.filtered.depth/filtered.F4.vcf.gz
OUTVCF2=~/030624_mergedVCF/Ionops_only_sorted_filtered.F4.vcf.gz
##############
##use bcftools sort and bcftools index to sort and index the vcf files to be merged, respectively
#bcftools sort $INVCF1 -Oz -o $OUTVCF1

#bcftools sort $INVCF2 -Oz -o $OUTVCF2

##now index
#bcftools index -t $OUTVCF1

#bcftools index -t $OUTVCF2
##########

###########
##use bcftools merge to merge the VCF files
#bcftools merge $OUTVCF1 $OUTVCF2 --force-samples -Oz -o ~/030624_mergedVCF/030624_Ionops_UKDipsTetsHex_merged.vcf.gz

#module unload bcftools-uoneasy/1.18-GCC-13.2.0

#################
##load GATK for SelectVariants (select biallelic loci only)
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

############
gatk SelectVariants -V $OUTVCF2 \
	--restrict-alleles-to BIALLELIC \
	--select-type-to-include SNP \
	-O ~/030624_mergedVCF/050624_IonopsOnly_BiSNP.vcf.gz
############

module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17


echo "DONE!!"

##lastline

