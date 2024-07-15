#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=16g
#SBATCH --time=01:00:00
#SBATCH --job-name=FixReadGroups
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile

##load Picard
module load picard-uoneasy/3.0.0-Java-17

################
##cd to the bam alignment files
cd ~/220524_alignments/bam_files/

meta=EKDL240001890-1A_222TKYLT4
#################
##fix read groups
# fix the read groups
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#    I=FLEET_2_${meta}.bam \
#    O=FLEET_2_${meta}_with_RG.bam \
#    SORT_ORDER=coordinate \
#    RGID=1A \
#    RGLB=EKDL24001890 \
#    RGPL=ILLUMINA \
#    RGPU=Unit1 \
#    RGSM=FLEET_2 \
#    CREATE_INDEX=True
###############
##fix read groups for HAM_1
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=HAM_1_${meta}.bam \
    O=HAM_1_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1B \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=Unit2 \
    RGSM=HAM_1 \
    CREATE_INDEX=True
##unload picard
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!"

##lastline
