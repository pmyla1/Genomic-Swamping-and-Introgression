#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=16g
#SBATCH --time=02:00:00
#SBATCH --job-name=select_dips_danica_anglica
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err

########
##source your profile
source $HOME/.bash_profile

#############
##load the GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##Index the vcf file using GATK IndexFeatureFile
gatk IndexFeatureFile -I /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_F4_133.ann.vcf.gz

##use GATK CreateSequenceDictionary to index our fasta file
gatk CreateSequenceDictionary -R /gpfs01/home/pmyla1/C_excelsa_V5.fasta -O /gpfs01/home/pmyla1/C_excelsa_V5.dict

##unload GATK module 
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
############

############
##load samtools module
module load samtools-uoneasy/1.18-GCC-12.3.0

##use samtools faidx to index the reference fasta file
samtools faidx /gpfs01/home/pmyla1/C_excelsa_V5.fasta -o /gpfs01/home/pmyla1/C_excelsa_V5.fasta.fai 

##unload samtools module
module unload samtools-uoneasy/1.18-GCC-12.3.0
############

#############
##load GATK module again
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

##use GATK SelectVariants to select only UK diploids with at least 2 individuals and all C. danica and C anglica
##restrict allleles to biallelic SNPs for Dsuite purposes with AF >0.50 to exclude rare variants
gatk SelectVariants -R /gpfs01/home/pmyla1/C_excelsa_V5.fasta -V /gpfs01/home/pmyla1/bcftools_gatk_output/Chrom1_F4_133.ann.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 --select "AF > 0.5 && AF < 1.0" \
 -sn CHA-01dl -sn CHA-02dl -sn JOR-01dl -sn JOR-02dl -sn JOR-03dl -sn JOR-04dl \
 -sn LAB-01dl -sn LAB-02dl -sn LAB-03dl -sn LAB-04dl -sn LAB-05dl -sn LAB-06dl -sn LAB-07dl -sn LAB-08dl \
 -sn NEN-01dl -sn NEN-02dl -sn NEN-03dl -sn NEN-04dl -sn NEN-05sl -sn NEN-06dl -sn NEN-07dl \
 -sn ODN-01dl -sn ODN-02dl -sn ODN-03dl -sn ODN-04dl -sn ODN-05dl -sn ODN-06dl -sn ODN-07dl \
 -sn JON-01hl -sn FOR-01hl -sn DAR-01hl -sn DAR-02hl -sn CUM-01hl -sn BRE-01hl -sn RYE-01hl -sn SCO-01hl \
 -sn SKF-01hl -sn SKF-02hl -sn SKF-03hl -sn SKF-04hl -sn SPU-01hl -sn SPU-02hl -sn SPU-03hl -sn SPU-04hl \
 -sn TET-01hl -sn TET-02hl -sn TET-03hl -sn TET-04hl -sn FRE-01hl \
 -O /gpfs01/home/pmyla1/bcftools_gatk_output/UK_dips_danica_all_bi_SNP_AF0.5.vcf.gz 

##unload the GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
###########

echo "DONE!!"

##lastline
