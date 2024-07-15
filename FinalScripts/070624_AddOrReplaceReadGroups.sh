#!/bin/bash

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=40g
#SBATCH --time=04:00:00
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
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=FLEET_2_${meta}.bam \
    O=FLEET_2_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1A \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=FLEET_2 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=HAM_1_${meta}.bam \
    O=HAM_1_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1B \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit2 \
    RGSM=HAM_1 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=LWS_${meta}.bam \
    O=LWS_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1C \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit3 \
    RGSM=LWS \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=NOT_${meta}.bam \
    O=NOT_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1D \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit4 \
    RGSM=NOT \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=PAR_2_${meta}.bam \
    O=PAR_2_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1E \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit5 \
    RGSM=PAR_2 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=Pen_1_${meta}.bam \
    O=Pen_1_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1F \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit6 \
    RGSM=Pen_1 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=SBAY_1_${meta}.bam \
    O=SBAY_1_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1G \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit7 \
    RGSM=SBAY_1 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=SPEY_2_${meta}.bam \
    O=SPEY_2_${meta}_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1H \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit8 \
    RGSM=SPEY_2 \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=Iac.bam \
    O=Iac_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1I \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit9 \
    RGSM=Iac \
    CREATE_INDEX=True
###############
# fix the read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I=Ime.bam \
    O=Ime_with_RG.bam \
    SORT_ORDER=coordinate \
    RGID=1J \
    RGLB=EKDL24001890 \
    RGPL=ILLUMINA \
    RGPU=unit10 \
    RGSM=Ime \
    CREATE_INDEX=True
###############

##unload picard
module unload picard-uoneasy/3.0.0-Java-17

echo "DONE!"

##lastline
