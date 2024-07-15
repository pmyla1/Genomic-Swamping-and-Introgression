#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=20g
#SBATCH --time=01:00:00
#SBATCH --job-name=Dquartets
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

##make an environmental variable for the VCF you want to use (050624_ld_pruned_allUkhex_someUKdips_someUKtets.vcf.gz)
VCF=~/050624_VCF/050624_ld_pruned_allUKhex_someUKdips_someUKtets.vcf.gz


##the trios with significantly elevated D-statistics and f4 admixture ratios
##GEO_BRE_SKF, SPU_GEO_BRE, FRE_ALO_BRE, FRE_GEO_BRE, FRE_LNL_BRE, TET_GEO_BRE, CUM_ALO_FRE, CUM_GEO_FRE, CUM_LNL_FRE, CUM_GEO_TET, JON_GEO_SPU,
##JON_GEO_TET, JON_GEO_SKF, JON_GEO_FRE, JON_AAH_FRE, JON_LNL_FRE, RYE_ALO_FRE, RYE_GEO_FRE

####use Dquartets (does NOT specify an outgroup) - use a jackknife block size of 2500
#Dsuite Dquartets -k 2500 -n 060624_JACKKNIFE_2500 $VCF ./Dquartets_SETS_SPECIES.txt

##RUN DQUARTETS AGAIN USING A DIFFERENT JACKKNIFE BLOCK SIZE TO SEE IF THE RESULTS ARE REPRODUCIBLE
##Jackknife block size of 5000
#Dsuite Dquartets -k 5000 -n 060624_JACKKNIFE_5000 $VCF ./Dquartets_SETS_SPECIES.txt

##Jackknife block size of 1000
#Dsuite Dquartets -k 1000 -n 060624_JACKKNIFE_1000 $VCF ./Dquartets_SETS_SPECIES.txt

##Jackknife block size of 10000
Dsuite Dquartets -k 10000 -n 060624_JACKKNIFE_10000 $VCF ./Dquartets_SETS_SPECIES.txt
#execute Dinvestigate on these trios 
#Dsuite Dinvestigate -w 50,25 -n 310524_Trios ~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz SIMULATED_SETs.txt test_trios.txt

module unload dsuite-uon/gcc11.3.0/0.5r57

echo "DONE!!"

##lastline
