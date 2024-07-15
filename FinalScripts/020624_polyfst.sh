#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=12g
#SBATCH --time=01:00:00
#SBATCH --job-name=sele_all_danica_anglica_UK_dips
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

###########
module load gcc-uoneasy/13.2.0 

cd ~/020624_polyfst_output/
########
##compile poly_fst.c script
#gcc ~/scripts/poly_fst.c -o poly_fst -lm
########
VCF=~/bcftools_gatk_output/ld_pruned_130524_allUKdips_allUKtets_allUKhex_copy.vcf

############
##GEO-hexaploid population contrast
##first for the GEO-SPU contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/SPU_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_SPU_polyfst.fst

##first for the GEO-TET contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/TET_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_TET_polyfst.fst

##first for the GEO-FRE contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/FRE_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_FRE_polyfst.fst

##first for the GEO-SKF contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/SKF_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_SKF_polyfst.fst

##first for the GEO-BRE contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/BRE_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_BRE_polyfst.fst

##first for the GEO-CUM contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/CUM_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_CUM_polyfst.fst

##first for the GEO-DAR contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/DAR_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_DAR_polyfst.fst

##first for the GEO-RYE contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/RYE_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_RYE_polyfst.fst

##first for the GEO-JON contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/JON_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_JON_polyfst.fst

##first for the GEO-SCO contrast
#./poly_fst -vcf $VCF -pop1 Populations/GEO_population.txt -pop2 Populations/SCO_population.txt -mis 0.8 -stat fst -out 0 > ./GEO_SCO_polyfst.fst
##############

################
##ALO-hexaploid pairwise population contrast
##first for the ALO-SKF contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/SCO_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_SKF_polyfst.fst
Â
##first for the ALO-SPU contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/SPU_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_SPU_polyfst.fst

##first for the ALO-TET contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/TET_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_TET_polyfst.fst

##first for the ALO-FRE contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/FRE_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_FRE_polyfst.fst

##first for the ALO-BRE contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/BRE_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_BRE_polyfst.fst

##first for the ALO-CUM contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/CUM_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_CUM_polyfst.fst

##first for the ALO-DAR contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/DAR_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_DAR_polyfst.fst

##first for the ALO-FOR contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/FOR_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_FOR_polyfst.fst

##first for the ALO-JON contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/JON_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_JON_polyfst.fst

##first for the ALO-RYE contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/RYE_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_RYE_polyfst.fst

##first for the ALO-SCO contrast
./poly_fst -vcf $VCF -pop1 Populations/ALO_population.txt -pop2 Populations/SCO_population.txt -mis 0.8 -stat fst -out 0 > ./ALO_SCO_polyfst.fst
#################

##################
##AAH-hexaploid population contrast
##first for the AAH-SKF contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/SKF_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_SKF_polyfst.fst

##first for the AAH-SPU contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/SPU_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_SPU_polyfst.fst

##first for the AAH-FRE contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/FRE_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_FRE_polyfst.fst

##first for the AAH-TET contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/TET_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_TET_polyfst.fst

##first for the AAH-CUM contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/CUM_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_CUM_polyfst.fst

##first for the AAH-BRE contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/BRE_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_BRE_polyfst.fst

##first for the AAH-DAR contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/DAR_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_DAR_polyfst.fst

##first for the AAH-FOR contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/FOR_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_FOR_polyfst.fst

##first for the AAH-JON contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/JON_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_JON_polyfst.fst

##first for the AAH-RYE contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/RYE_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_RYE_polyfst.fst

##first for the AAH-SCO contrast
./poly_fst -vcf $VCF -pop1 Populations/AAH_population.txt -pop2 Populations/SCO_population.txt -mis 0.8 -stat fst -out 0 > ./AAH_SCO_polyfst.fst
############

################
##LNL-hexaploid population contrast
##first for the LNL-SKF contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/SKF_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_SKF_polyfst.fst

##first for the LNL-SPU contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/SPU_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_SPU_polyfst.fst

##first for the LNL-FRE contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/FRE_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_FRE_polyfst.fst

##first for the LNL-TET contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/TET_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_TET_polyfst.fst

##first for the LNL-CUM contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/CUM_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_CUM_polyfst.fst

##first for the LNL-BRE contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/BRE_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_BRE_polyfst.fst

##first for the LNL-DAR contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/DAR_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_DAR_polyfst.fst

##first for the LNL-FOR contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/FOR_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_FOR_polyfst.fst

##first for the LNL-JON contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/JON_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_JON_polyfst.fst

##first for the LNL-RYE contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/RYE_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_RYE_polyfst.fst

##first for the LNL-SCO contrast
./poly_fst -vcf $VCF -pop1 Populations/LNL_population.txt -pop2 Populations/SCO_population.txt -mis 0.8 -stat fst -out 0 > ./LNL_SCO_polyfst.fst


module unload gcc-uoneasy/13.2.0 

echo "DONE!!"

##lastline
