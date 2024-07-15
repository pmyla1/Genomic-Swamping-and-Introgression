#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=12g
#SBATCH --time=01:00:00
#SBATCH --job-name=bcftools_query
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err

###################
##source bash profile
source $HOME/.bash_profile

##Use bcftools to extract only the coordinates found in the C_excelsa_V5.fasta reference
##load bcftools module 
module load bcftools-uoneasy/1.18-GCC-13.2.0

###
##check the sample names in the vcf file with bcftools query -l 
bcftools query -l /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz > /gpfs01/home/pmyla1/bcftools_gatk_output/100524_sample_names.txt

####
##extract specific columns from the vcf file
bcftools query -f '%CHROM %POS %REF %ALT\n' /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz -o /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom_pos_ref_alt_100524.txt

####
##extract only chromosome 1 (Cexcelsa_scaf_1) because this is the only scaffold in the C_excelsa_V5.fasta reference
#bcftools view /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz --regions Cexcelsa_scaf_1 -o /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_F4_133.ann.vcf.gz

##extract the chromosome 1 coordinates in the C_excelsa_V5.fasta reference 
bcftools view /gpfs01/home/pmyla1/reheadered.F4_133.ann.vcf.gz --regions Cexcelsa_scaf_1:1-1407794 -o /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_coords_1_1407794_F4_133.ann.vcf.gz 

##now extract this information for UK diploids, UK tetraploids, all C. danica, and all putative C. anglica
##use grep to extract the samples you want from 100524_samples_names.txt
grep 'AAH\|ROT\|SKN\|ALO\|ELI\|ERS\|LAL\|FTW\|LNL\|LOS\|NEI\|SCU\|GEO\|BNK\|CHA\|JOR\|LAB\|NEN\|ODN\|JON\|FOR\|DAR\|CUM\|BRE\|RYE\|SCO\|SKF\|SPU\|TET\|FRE' /gpfs01/home/pmyla1/bcftools_gatk_output/100524_sample_names.txt > /gpfs01/home/pmyla1/bcftools_gatk_output/UK_dips_tets_danica_anglica_samples.txt

bcftools query -f '%CHROM %POS %REF %ALT\n' /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_coords_1_1407794_F4_133.ann.vcf.gz -S /gpfs01/home/pmyla1/bcftools_gatk_output/UK_dips_tets_danica_anglica_samples.txt -o /gpfs01/home/pmyla1/bcftools_gatk_output/UK_dips_tets_danica_anglica_samples.vcf.gz 
##unload module
module unload bcftools-uoneasy/1.18-GCC-13.2.0
###############

#################
##load GATK module
#module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##index the new VCF file for the co-ordinates chosen above (Cexcelsa_scaf1:1-1407794)
#gatk IndexFeatureFile -F /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_coords_1_1407794_F4_133.ann.vcf.gz

##unload GATK module again
#module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
###############

###############
##Indexing C_excelsa_V5.fasta
##load samtools module for indexing
#module load samtools-uoneasy/1.18-GCC-12.3.0

##use samtools faidx to index the reference fasta file
#samtools faidx /gpfs01/home/pmyla1/C_excelsa_V5.fasta -o /gpfs01/home/pmyla1/C_excelsa_V5.fasta.fai 

##unload samtools module
#module unload samtools-uoneasy/1.18-GCC-12.3.0
############

#############
##Selecting the UK diploids, all C. danica, and all putatitve C. anglica
##load GATK module again
#module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##use GATK SelectVariants to select only UK diploids with at least 2 individuals and all C. danica and C anglica
##restrict allleles to biallelic SNPs for Dsuite purposes with AF >0.50 to exclude rare variants
#gatk SelectVariants -R /gpfs01/home/pmyla1/C_excelsa_V5.fasta -V /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_F4_133.ann.vcf.gz \
# --select-type-to-include SNP \
# --restrict-alleles-to BIALLELIC \
# --select "AF > 0.5 && AF < 1.0" \
# -sn CHA-01dl -sn CHA-02dl -sn JOR-01dl -sn JOR-02dl -sn JOR-03dl -sn JOR-04dl \
# -sn LAB-01dl -sn LAB-02dl -sn LAB-03dl -sn LAB-04dl -sn LAB-05dl -sn LAB-06dl -sn LAB-07dl -sn LAB-08dl \
# -sn NEN-01dl -sn NEN-02dl -sn NEN-03dl -sn NEN-04dl -sn NEN-05sl -sn NEN-06dl -sn NEN-07dl \
# -sn ODN-01dl -sn ODN-02dl -sn ODN-03dl -sn ODN-04dl -sn ODN-05dl -sn ODN-06dl -sn ODN-07dl \
# -sn JON-01hl -sn FOR-01hl -sn DAR-01hl -sn DAR-02hl -sn CUM-01hl -sn BRE-01hl -sn RYE-01hl -sn SCO-01hl \
# -sn SKF-01hl -sn SKF-02hl -sn SKF-03hl -sn SKF-04hl -sn SPU-01hl -sn SPU-02hl -sn SPU-03hl -sn SPU-04hl \
# -sn TET-01hl -sn TET-02hl -sn TET-03hl -sn TET-04hl -sn FRE-01hl \
# -O /gpfs01/home/pmyla1/bcftools_gatk_output/UK_dips_danica_all_bi_SNP_AF0.5.vcf.gz 

##unload GATK module
#module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##############

echo "DONE!!"

##lastline
