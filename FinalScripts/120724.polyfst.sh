#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=06:00:00
#SBATCH --job-name=polyfst
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load gcc 
module load gcc-uoneasy/13.2.0  

##cd to appropriate directory
cd ~/120724.polyfst
##poly_fst calculating fst
#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/officinalis.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/pyrenaica.officinalis.fst.txt

#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/anglica.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/pyrenaica.anglica.fst.txt

#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/danica.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/pyrenaica.danica.fst.txt

#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/officinalis.txt -pop2 ~/120724.polyfst/danica.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/officinalis.danica.fst.txt

#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/officinalis.txt -pop2 ~/120724.polyfst/anglica.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/officinalis.anglica.fst.txt

#~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/anglica.txt -pop2 ~/120724.polyfst/danica.txt -mis 0.9 -maf 0.05 -stat fst > ~/120724.polyfst/anglica.danica.fst.txt
################
##polyfst calculating dxy = nucleotide diversity
~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/officinalis.txt -mis 0.9 -maf 0.05 -stat dxy > ~/pyrenaica.officinalis.dxy

~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/anglica.txt -mis 0.9 -maf 0.05 -stat dxy > ~/pyrenaica.anglica.dxy

~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/pyrenaica.txt -pop2 ~/120724.polyfst/danica.txt -mis 0.9 -maf 0.05 -stat dxy > ~/pyrenaica.danic.dxy

~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/danica.txt -pop2 ~/120724.polyfst/officinalis.txt -mis 0.9 -maf 0.05 -stat dxy > ~/danica.officinalis.dxy

~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/anglica.txt -pop2 ~/120724.polyfst/officinalis.txt -mis 0.9 -maf 0.05 -stat dxy > ~/anglica.officinalis.dxy

~/scripts/poly_fst -vcf ~/120624_LD.Pruned.Ionops.allUKsamples.copy.vcf -pop1 ~/120724.polyfst/anglica.txt -pop2 ~/120724.polyfst/danica.txt -mis 0.9 -maf 0.05 -stat dxy > ~/anglica.danica.dxy
##############

##unload module
module unload gcc-uoneasy/13.2.0  

echo "DONE!!"

##lastline
