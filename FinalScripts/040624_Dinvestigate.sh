#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20g
#SBATCH --time=01:00:00
#SBATCH --job-name=Dsuite
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
#######
##This script was written by Luke Archer (2024) and can be used to perform Dsuite for a variety of different population
##combinations, including UK diploids, UK tetraploids, UK hexaploids and C. danica

source $HOME/.bash_profile

##load Dsuite module
module load dsuite-uon/gcc11.3.0/0.5r57

cd ~/Dsuite_files/

##make an environmental variable for the VCF you want to use (ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz)
#VCF=~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz


##the trios with significantly elevated D-statistics and f4 admixture ratios
##GEO_BRE_SKF, SPU_GEO_BRE, FRE_ALO_BRE, FRE_GEO_BRE, FRE_LNL_BRE, TET_GEO_BRE, CUM_ALO_FRE, CUM_GEO_FRE, CUM_LNL_FRE, CUM_GEO_TET, JON_GEO_SPU,
##JON_GEO_TET, JON_GEO_SKF, JON_GEO_FRE, JON_AAH_FRE, JON_LNL_FRE, RYE_ALO_FRE, RYE_GEO_FRE

#execute Dinvestigate on these trios with a window size of 10 SNPs and 1 SNP window step size for higher resolution 
#Dsuite Dinvestigate -w 10,1 -n 310524_Trios ~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz SIMULATED_SETS.txt test_trios.txt
Dsuite Dinvestigate -w 50,1 -n 310524_Trios ~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz SIMULATED_SETS.txt test_trios.txt


module unload dsuite-uon/gcc11.3.0/0.5r57

echo "DONE!!"

##lastline
