#!/bin/bash

#source bash profile
source $HOME/.bash_profile

##utilise Dsuite
cd ~/Dsuite/

##execute Dtrios with a jackknife number of 50 (splits the SNPs into 50 blocks)
./Build/Dsuite Dtrios -k 50 -t ../Desktop/Individual_Project_files/Dsuite_files/150524_experimental_TREE_FILE.nwk --ABBAclustering -o 150524_EU_UK_ld_pruned ../Desktop/Individual_Project_files/VCF_files/150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf.gz ../Desktop/Individual_project_files/Dsuite_files/150524_SETs.txt


##now execute Fbranch
./Build/Dsuite Dinvestigate -w 100,50 ../Desktop/Individual_Project_files/VCF_files/150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf.gz ../Desktop/Individual_project_files/Dsuite_files/150524_SETs.txt ../Desktop/Individual_project_files/Dsuite_files/test_trios.txt

##now execute Fbranch
./Build/Dsuite Fbranch -Z ../Desktop/Individual_project_files/Dsuite_files/TREE_FILE.nwk ../Desktop/Individual_project_files/Dsuite_files/150524_EU_UK_ld_pruned_tree.txt

./utils/dtools.py --outgroup Outgroup C_danica_C_anglica_C_officinalis_localFstats__100_50.txt ../Desktop/Individual_project_files/Dsuite_files/TREE_FILE.nwk


#######160524 perform Dquartets (NO OUTGROUP)
./Build/Dsuite Dquartets -k 50 -t ../Desktop/Individual_Project_files/Dsuite_files/160524_Dquartet_TREE_FILE.nwk ../Desktop/Individual_Project_files/VCF_files/150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf.gz ../Desktop/Individual_project_files/Dsuite_files/160524_Dquartet_SETs.txt
