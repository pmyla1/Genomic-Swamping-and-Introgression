#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=select_UK_hexaploids_only
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##module load gatk
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

###########
##use GATK to select only the UK hexaploids from the original VCF file
gatk SelectVariants \
 -V ~/reheadered.F4_133.ann.vcf.gz \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 -sn BRE_1 -sn CUM_1 -sn DAR_1 -sn DAR_3 -sn FOR_1 -sn FRE_013 \
 -sn JON_001 -sn RYE_1 -sn SCO_1 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 \
 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 -O ~/140524_UKhexaploids_only_vcf/140524_UKhexaploids_only.vcf.gz

##module unload
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17
##########

############
module load samtools-uoneasy/1.18-GCC-12.3.0
##make a copy of the vcf and then unzip 
cp ~/140524_UKhexaploids_only_vcf/140524_UKhexaploids_only.vcf.gz ~/140524_UKhexaploids_only_vcf/140524_UKhexaploids_only_copy.vcf.gz

gunzip ~/140524_UKhexaploids_only_vcf/140524_UKhexaploids_only_copy.vcf.gz

module unload samtools-uoneasy/1.18-GCC-12.3.0
##########

##########
##now ld prune the VCF file using prune_ld.c
module load gcc-uoneasy/13.2.0
##configure prune_ld
gcc ~/scripts/prune_ld.c -o ~/140524_UKhexaploids_only_vcf/140524_prune_ld -lm

#########
cd ~/140524_UKhexaploids_only_vcf/
##execute 140524_prune_ld on the 140524_UKhexaploids_only.vcf.gz
./140524_prune_ld -vcf ./140524_UKhexaploids_only_copy.vcf -mis 0.8 -maf 0.05 -r2 100 50 0.1 > ./140524_ld_pruned_UKhexaploids_only_copy.vcf 

module unload gcc-uoneasy/13.2.0
##########

#########
##bgzip the 140524_ld_pruned_UKhexaploids_only_copy.vcf 
cp ./140524_ld_pruned_UKhexaploids_only_copy.vcf ./140524_ld_pruned_UKhexaploids_only.vcf

module load htslib-uoneasy/1.18-GCC-13.2.0

##bgzip the ld_pruned.vcf file
bgzip ./140524_ld_pruned_UKhexaploids_only.vcf

module unload htslib-uoneasy/1.18-GCC-13.2.0
##########

echo "DONE!!"

##lastline
