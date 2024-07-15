#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=02:00:00
#SBATCH --job-name=select_chrom1_coords
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load gcc module
module load gcc-uoneasy/13.2.0

##compile poly_fst script
gcc /gpfs01/home/pmyla1/scripts/poly_fst.c -o /gpfs01/home/pmyla1/bcftools_gatk_output/poly_fst -lm

###########
##change directory to bcftools_gatk_output
cd /gpfs01/home/pmyla1/bcftools_gatk_output/

##use grep to pull the hexaploid C. danica and C. anglica populations you want to contrast from the 140524_inds.txt file 
grep "SKF" ./140524_inds.txt > ./SKF_C_anglica.txt

grep "SPU" ./140524_inds.txt > ./SPU_C_anglica.txt

grep "TET" ./140524_inds.txt > ./TET_C_anglica.txt

grep "FRE" ./140524_inds.txt > ./FRE_C_anglica.txt

grep "BRE" ./140524_inds.txt > ./BRE_C_danica.txt

grep "CUM" ./140524_inds.txt > ./CUM_C_danica.txt

grep "DAR" ./140524_inds.txt > ./DAR_C_danica.txt

grep "FOR" ./140524_inds.txt > ./FOR_C_danica.txt

grep "JON" ./140524_inds.txt > ./JON_C_danica.txt

grep "RYE" ./140524_inds.txt > ./RYE_C_danica.txt

grep "SCO" ./140524_inds.txt > ./SCO_C_danica.txt
#############

##############
##unzip the ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz
cp ./ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf.gz

gunzip ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf.gz

##now use poly_fst on different hexaploid populations comparing the C_anglica with the C_danica populations
mkdir SKF_poly_fst_output
##first for the SKF vs BRE contrast
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./BRE_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_BRE_C_dan_C_ang_constrast.fst

##now for SKF vs CUM
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./CUM_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_CUM_C_dan_C_ang_constrast.fst

##now for SKF vs DAR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./DAR_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_DAR_C_dan_C_ang_constrast.fst

##now for SKF vs FOR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./FOR_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_FOR_C_dan_C_ang_constrast.fst

##now for SKF vs JON
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./JON_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_JON_C_dan_C_ang_constrast.fst

##now for SKF vs RYE 
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./RYE_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_RYE_C_dan_C_ang_constrast.fst

##now for SKF vs SCO
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SKF_C_anglica.txt -pop2 ./SCO_C_danica.txt -mis 0.9 -stat fst > ./SKF_poly_fst_output/SKF_SCO_C_dan_C_ang_constrast.fst
##################

##################
mkdir SPU_poly_fst_output
##now for SPU vs BRE
##first for the SKF vs BRE contrast
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./BRE_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_BRE_C_dan_C_ang_constrast.fst

##now for SPU vs CUM
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./CUM_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_CUM_C_dan_C_ang_constrast.fst

##now for SPU vs DAR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./DAR_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_DAR_C_dan_C_ang_constrast.fst

##now for SPU vs FOR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./FOR_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_FOR_C_dan_C_ang_constrast.fst

##now for SPU vs JON
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./JON_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_JON_C_dan_C_ang_constrast.fst

##now for SPU vs RYE 
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./RYE_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_RYE_C_dan_C_ang_constrast.fst

##now for SPU vs SCO
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./SPU_C_anglica.txt -pop2 ./SCO_C_danica.txt -mis 0.9 -stat fst > ./SPU_poly_fst_output/SPU_SCO_C_dan_C_ang_constrast.fst
#################

################
mkdir TET_poly_fst_output
##now for TET vs BRE contrast
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./BRE_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_BRE_C_dan_C_ang_constrast.fst

##now for TET vs CUM
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./CUM_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_CUM_C_dan_C_ang_constrast.fst

##now for TET vs DAR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./DAR_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_DAR_C_dan_C_ang_constrast.fst

##now for TET vs FOR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./FOR_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_FOR_C_dan_C_ang_constrast.fst

##now for TET vs JON
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./JON_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_JON_C_dan_C_ang_constrast.fst

##now for TET vs RYE 
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./RYE_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_RYE_C_dan_C_ang_constrast.fst

##now for TET vs SCO
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./TET_C_anglica.txt -pop2 ./SCO_C_danica.txt -mis 0.9 -stat fst > ./TET_poly_fst_output/TET_SCO_C_dan_C_ang_constrast.fst
#################

################
mkdir FRE_poly_fst_output
##now for FRE vs BRE contrast
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./BRE_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_BRE_C_dan_C_ang_constrast.fst

##now for FRE vs CUM
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./CUM_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_CUM_C_dan_C_ang_constrast.fst

##now for FRE vs DAR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./DAR_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_DAR_C_dan_C_ang_constrast.fst

##now for FRE vs FOR
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./FOR_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_FOR_C_dan_C_ang_constrast.fst

##now for FRE vs JON
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./JON_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_JON_C_dan_C_ang_constrast.fst

##now for FRE vs RYE 
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./RYE_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_RYE_C_dan_C_ang_constrast.fst

##now for FRE vs SCO
./poly_fst -vcf ./ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf -pop1 ./FRE_C_anglica.txt -pop2 ./SCO_C_danica.txt -mis 0.9 -stat fst > ./FRE_poly_fst_output/FRE_SCO_C_dan_C_ang_constrast.fst
###############

##unload gcc module
module unload gcc-uoneasy/13.2.0

echo "DONE!!!"

##lastline
