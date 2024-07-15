#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=24g
#SBATCH --time=06:00:00
#SBATCH --job-name=FastaAlternateReferenceMaker
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##source bash profile
source $HOME/.bash_profile

##load GATK module
module load gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

######
##environmental variables for the reference (REF), output (OUT), and vcf (VCF)
REF=~/280524_C_excelsa_V5_reference/C_excelsa_V5.fasta
#OUTJAR1=~/110724.GenesOfInterest/110724.g50778.JAR1.coords.fasta
#OUTAVI2=~/110724.GenesOfInterest/110724.g50851.AVI2.coords.fasta
OUTGRP7=~/110724.GenesOfInterest/110724.g50878.GRP7.coords.fasta
VCF=~/120624_LD.Pruned.Ionops.allUKsamples.vcf.gz
######

##use GATK IndexFeatureFile to index the VCF file
gatk IndexFeatureFile -I $VCF

##use GATK FastaAlternateReferenceMaker to extract the CDS from the g50778.t1 gene encoding JAR1 (Jasmonate Resistant 1)
# gatk FastaAlternateReferenceMaker \
 #  -R $REF \
  # -O $OUTJAR1 \
   #-L Cexcelsa_scaf_6:5890772-5890782 -L Cexcelsa_scaf_6:5891040-5891069 -L Cexcelsa_scaf_6:5891770-5892090 \
  # -L Cexcelsa_scaf_6:5892184-5892285 -L Cexcelsa_scaf_6:5892388-5893175 -L Cexcelsa_scaf_6:5893260-5893801 \
  # -V $VCF \
###########
##use GATK FastaAlternateReferenceMaker to extract the CDS from the g50851.t1 gene encoding AVI2
 #gatk FastaAlternateReferenceMaker \
 #  -R $REF \
 #  -O $OUTAVI2 \
 #  -L Cexcelsa_scaf_6:6172087-6172589 -L Cexcelsa_scaf_6:6173101-6173130 -L Cexcelsa_scaf_6:6173309-6173426 \
 #  -L Cexcelsa_scaf_6:6173511-6173603 -L Cexcelsa_scaf_6:6173925-6174062 -L Cexcelsa_scaf_6:6174172-6174299 \
 #  -L Cexcelsa_scaf_6:6174397-6174557 -L Cexcelsa_scaf_6:6174657-6174747 -L Cexcelsa_scaf_6:6174841-6174897 \
 #  -L Cexcelsa_scaf_6:6174972-6175092 \
 #  -V $VCF \
###########
gatk FastaAlternateReferenceMaker \
	-R $REF \
	-O $OUTGRP7 \
	-L Cexcelsa_scaf_6:6271476-6271865 -L Cexcelsa_scaf_6:6272148-6272261 \
	-V $VCF
##unload GATK module
module unload gatk-uoneasy/4.4.0.0-GCCcore-12.3.0-Java-17

echo "DONE!!"

##lastline

