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

cd ~/230624_Dsuite/

##make an environmental variable for the VCF you want to use (120624_LD.pruned.allUKsamples.vcf.gz)
VCF=~/120624_LD.Pruned.Ionops.allUKsamples.vcf.gz

###########
##execute Dquartets on the dataset to obtain the assumed relationships between the species excluding the Ionopsidium outgroup samples
#Dsuite Dquartets -k 4000 -o 120624_experimental $VCF SETS_SPECIES.txt

#execute Dtrios using the appropriate outgroup (Iac, Ime, Iab_1, Iab_2) using 4000 Jack-knife blocks (-k 4000)
#without a tree
#Dsuite Dtrios -k 4000 -o 240624_Dtrios --ABBAclustering $VCF SETS_SPECIES.txt

#############
##now use Dinvestigate to investigate the trios with elevated D-statistics
##first for the anglica officinalis     danica trio
#Dsuite Dinvestigate -w 50,25 -n 50_25 $VCF SETS_SPECIES.txt 240624_testtrios.txt

##now use different SNP window sizes and step sizes 
#Dsuite Dinvestigate -w 100,10 -n 100_10_pyr_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
## try for 100SNP windows and 25 SNP step size
#Dsuite Dinvestigate -w 100,25 -n 100_25_pyr_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
##try for 50SNP windows and 25 SNP step size
#Dsuite Dinvestigate -w 50,25 -n 50_25_pyr_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
###############

#############
##first for the anglica	officinalis	danica trio
#Dsuite Dinvestigate -w 50,25 -n 50_25_ang_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt

#Dsuite Dinvestigate -w 100,10 -n 100_10_ang_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt

#Dsuite Dinvestigate -w 100,25 -n 100_25_ang_off_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
##########
#next for the pyrenica	anglica	danica trio
#Dsuite Dinvestigate -w 100,25 -n 100_25_pyr_ang_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt

Dsuite Dinvestigate -w 500,25 -n 100_25_pyr_ang_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt

Dsuite Dinvestigate -w 1000,25 -n 100_25_pyr_ang_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
#Dsuite Dinvestigate -w 100,10 -n 100_10_pyr_ang_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt

#Dsuite Dinvestigate -w 50,25 -n 50_25_pyr_ang_dan $VCF SETS_SPECIES.txt 240624_testtrios.txt
###########

##################
##now for the pyrenica	officinalis anglica trio
#Dsuite Dinvestigate -w 100,10 -n 100_10_pyr_off_ang $VCF SETS_SPECIES.txt 240624_testtrios.txt

#Dsuite Dinvestigate -w 50,25 -n 50_25_pyr_off_ang $VCF SETS_SPECIES.txt 240624_testtrios.txt

#Dsuite Dinvestigate -w 100,25 -n 100_25_pyr_off_ang $VCF SETS_SPECIES.txt 240624_testtrios.txt
##################

###############
##now execute Fbranch on the _tree.txt output from Dtrios
#Dsuite Fbranch -Z 120624_Tree.nwk 120624_Dtrios_tree.txt > ~/110624_Dsuite/120624_fbranch.txt
########
##set up the dtools.py script for plotting fbranch

ø
module unload dsuite-uon/gcc11.3.0/0.5r57

echo "DONE!!"

##lastline
