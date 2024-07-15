#!/bin/bash

#SBATCH --partition=shortq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=12g
#SBATCH --time=01:00:00
#SBATCH --job-name=MAPREADSSAMTOOLSSORT
#SBATCH --output=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.out
#SBATCH --error=/gpfs01/home/pmyla1/slurm_output_error/slurm-%x-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmyla1@exmail.nottingham.ac.uk

##################
##setup 
source $HOME/.bash_profile
##CONVERT SAM FILES TO BAM FILES, SORT BAMS, INDEX, THEN FLAGSTAT WITH SAMTOOLS
module load samtools-uoneasy/1.18-GCC-12.3.0

##change directory to the 220524_alignments folder
cd ~/220524_alignments/

mkdir -p bam_files/

metadata=EKDL240001890-1A_222TKYLT4
##########
##samtools view with 4 threads, output headers(-h), in bam format (-b)
##samtools sort with 4 threads, using .bam as input
##samtool index
##samtools flagstat on the sorted bam
################
#on FLE_2
#samtools view -@ 4 -h -b ./FLEET_2_${metadata}_aln-pe.sam -o ./bam_files/FLEET_2_${metadata}.bam
#samtools sort -@ 4 -o ./bam_files/FLEET_2_${metadata}.sorted.bam ./bam_files/FLEET_2_${metadata}.bam
#samtools index ./bam_files/FLEET_2_${metadata}.sorted.bam
#samtools flagstat ./bam_files/FLEET_2_${metadata}.sorted.bam > ./bam_files/FLEET_2_${metadata}.flagstats
##now on HAM_1
#samtools view -@ 4 -h -b ./HAM_1_${metadata}_aln-pe.sam -o ./bam_files/HAM_1_${metadata}.bam
#samtools sort -@ 4 -o ./bam_files/HAM_1_${metadata}.sorted.bam ./bam_files/HAM_1_${metadata}.bam
#samtools index ./bam_files/HAM_1_${metadata}.sorted.bam
#samtools flagstat ./bam_files/HAM_1_${metadata}.sorted.bam > ./bam_files/HAM_1_${metadata}.flagstats
#############
##now on LWS
samtools view -@ 4 -h -b ./LWS_${metadata}_aln-pe.sam -o ./bam_files/LWS_${metadata}.bam
samtools sort -@ 4 -o ./bam_files/LWS_${metadata}.sorted.bam ./bam_files/LWS_${metadata}.bam
samtools index ./bam_files/LWS_${metadata}.sorted.bam
samtools flagstat ./bam_files/LWS_${metadata}.sorted.bam > ./bam_files/LWS_${metadata}.flagstats

module unload samtools-uoneasy/1.18-GCC-12.3.0

echo "DONE!!"

##lastline
