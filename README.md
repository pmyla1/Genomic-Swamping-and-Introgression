# Genomic-Swamping-and-Introgression
This repository should allow the user to reproduce an analysis of the extent of **introgression** between UK accessions of invasive hexaploid ***Cochlearia danica*** and native tetraploid species ***Cochlearia officinalis***, utilising tools such as **Dsuite** and **Twisst**. 

# Background

The ***Cochlearia* species complex** is a promising study system for the **evolutionary genomics of adaptation** due to the rapid acquisition of traits such as **cold tolerance**, **salt tolerance**, and **heavy-metal tolerance**, in a short evolutionary timescale ([Wolf et al., 2021](https://doi.org/10.7554/eLife.71572)). This species complex **thrives** in cold **Alpine and Arctic** environments which contrasts greatly from the preferred **Mediterranean habitat** of its sole **sister taxa *Ionopsidium*** (Wolf et al., 2021). The genus *Cochlearia* consists of **16 accepted species** and **4 subspecies**, and started to **rapidly diversify** due to periodic **climatic fluctuations** that occurred during the **Middle** and **Late Pleistocene** (between **0.77-0.012** million years ago; Wolf et al., 2021).

The genus *Cochlearia* displays a **wide range of cytotypes** and levels of ploidy, ranging from **diploids** to tetraploids to hexaploids to **octaploids**, which makes this species complex an interesting **model species** for the study of **adaptations to whole genome duplication** ([Bray et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.017939v1.full)). Since whole genome duplication and polyploidy are **major effect mutations**, many **cellular**, **ionomic**, and **molecular** processes are **disrupted**, especially those pertaining to **sister chromatid segregation** during meiosis, **DNA repair**, and **recombination**, amongst others ([Yant & Schmickl, 2021](https://pubmed.ncbi.nlm.nih.gov/33454987/)). Therefore, **neopolyploids** must overcome these initial challenges associated with sister chromatid segregation during meiosis, to enable **reduced crossover numbers** and to **prevent chromosomal breakages** during anaphase (Bray et al., 2020). Some of the **genes** found to be **under selection in neopolyploids** relative to their diploid counterparts are involved in biological processes such as **DNA repair**, **recombination**, **sister chromatid segregation**, amongst others, however, despite this **process-level convergence** there appears to be **low orthologue-level convergenc** ([Bray et al., 2023](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v1.full)). 

*Cochlearia danica* is a **highly invasive, salt tolerant** species in the *Cochlearia* genus that is **native** to the **Atlantic coasts of Europe** in countries including Denmark, Belgium, the UK, Germany, Finland, Ireland, Norway, Russia, and Sweden amongst others (Fekete et al., 2018). Despite being native to the Atlantic coasts of Europe, *C. danica* has **spread rapidly throughout Central Europe** since the 1970s, especially **along roadsides** which is attributed to the **widespread use of de-icing salts** ([Fekete et al., 2018](http://dx.doi.org/10.23855/preslia.2018.023)). The **rate of spread** of *C. danica* along Central European roadsides was estimated to be **approximately 62-65km/year**, and the soil in which it grows and thrives is characterised by **high salt content** (Fekete et al., 2018). Interestingly, this species has been reported to undergo **rapid and marked/remarkable** changes in **population size**, and one Hungarian population was found to decrease in size by 99% between 2016 and 2017 (Fekete et al., 2018).

Introgression is the exchange of genetic material between species that results from hybridization and recurrent backcrossing ([Wang et al., 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10504873/)).

The aim of this project is to detect and quantify the extent of gene flow/introgression between the invasive hexaploid halophyte *C. danica* and native UK populations of *C. anglica* (6X) or UK populations of *C. officinalis* (4X) by calculating ABBA-BABA statistics on genome-wide SNP data using software such as Dsuite ([Malinsky, 2021](https://doi.org/10.1111/1755-0998.13265)). 


# Installation of Software and Dependencies

## Dsuite
[Dsuite](https://github.com/millanek/Dsuite) is a software program developed to quickly calculate Patterson's D (ABBA-BABA), and the f4-ratio statistics across many populations and/or species. 

This software takes a VCF file and an explicitly stated phylogenetic tree in the Newick format as input and uses "parsimony informative" single nucleotide polymorphisms from a quartet of species to detect the occurrence of gene flow between species on internal branches of the tree. 

Citation: Malinsky, M., Matschiner, M. and Svardal, H. (2021) Dsuite ‐ fast D‐statistics and related admixture evidence from VCF files. Molecular Ecology Resources 21, 584–595. doi:[https://doi.org/10.1111/1755-0998.13265](https://doi.org/10.1111/1755-0998.13265)

To install the program on macOS run the following commands for the main program:

```
git clone https://github.com/millanek/Dsuite.git
cd Dsuite
make
```
In order to execute Dsuite commands (e.g. Dtrios), you can navigate to the Build folder and run the Dsuite executable with the following command `./Build/Dsuite` which shows the available commands. To execute the Dtrios command you can type `./Build/Dsuite Dtrios`.

I executed the Dsuite Dtrios command locally utilising the following commands:
```
./Build/Dsuite Dtrios -t ../Desktop/Individual_Project_files/Dsuite_files/TREE_FILE.nwk -o 120524_Dsuite_stats --ABBAclustering ../Desktop/Individual_Project_files/ld_pruned_110524_WG_allUKhex_allUKdips_someUKtets_copy.vcf ../Desktop/Individual_Project_files/Dsuite_files/SETs.txt 
```

The TREE_FILE.nwk was structured as following, where the Outgroup was the UK diploid species *Cochlearia pyrenaica* because data for *Ionopsidium* (the true sister taxa and real outgroup) was not available, and UK diploids were indicated as the alternative outgroup by `Dsuite Dquartets`:
```
(Outgroup,(officinalis,(anglica,danica)));
```
The SETs.txt file has the following structure with the individual ID and the group ID (i.e. the species) separated by a tab, and is demonstrated below:
```
BNK_21 Outgroup
CHA_1  Outgroup
CHA_2  Outgroup
JOR_1  Outgroup
...
AAH_1  officinalis
AAH_2  officinalis
AAH_3  officinalis
AAH_4  officinalis
...
BRE_1  danica
CUM_1  danica
DAR_1  danica
DAR_3  danica
...
SKF_002  anglica
SKF_003  anglica
SKF_005  anglica
SKF_009  anglica
```

# Data Generation (17/05/2024)

Additional *C. danica* and *Ionopsidium* Illumina paired-end sequencing data was provided by Yant (2024), and the sequencing reads (in fq.gz/fastq.gz format) were processed loosely following the guidelines outlined in [ngs_pipe](https://github.com/mattheatley/ngs_pipe/blob/main/README) from Healey (2024).

## FastQC and MultiQC for sequencing quality control reports

An example command for producing a fastqc report for a single population (FLEET_2 in this example) can be found below:
```
###load fastqc for performing sequencing quality control reports
module load fastqc-uoneasy/0.12.1-Java-11

##change directory to the fastq.gz files 
cd /file/path/to/fastq.gz files

##perform fastqc on the FLEET_2 sequencing reads
fastqc -o ../170524_fastqc/ ./FLE_2/*.fq.gz

##REPEAT for the other populations
```

Subsequently, a MultiQC report can be generated by utilising the directory containing the FastQC reports generated in the previous stage as input files (i.e. the .fastqc.zip files). The command to produce the MultiQC report can be found below.
```
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a

##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc /gpfs01/home/pmyla1/170524_fastqc/.*fastqc.zip
```
The MultiQC plots and reports include a directory for png, svg, or pdf versions of the plots.

## 22/05/2024 - Whole Pipeline from Data Generation to VCF

The Illumina paired-end sequencing data provided by Yant (2024) were processed following the steps outlined in the [ngs_pipe](https://github.com/mattheatley/ngs_pipe/blob/main/README) README page written by Healey (2024). 

# Stage 1: Trimming the adapters from Illumina paired-end sequencing reads

Firstly, Nextera Transposase adapters were trimmed from the reads with Trimmomatic (version 0.39), by specifying `ILLUMINACLIP:NexteraPE-PE.fa:2:40:15` on the command line. Reads with PHRED scores < 20 and a minimum length of 25 were trimmed ((`SLIDINGWINDOW:4:20` & `MINLEN:25`, respectively). An example trimmomatic command for one population (HAM_1) can be found below:

```
##Trimmomatic - trims Nextera transposase adapters from the Illumina reads
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
        PE -phred33 ../HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ../HAM_1/HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz \
        ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.trimmed.fq.gz ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_1.orhpan.fq.gz \
        ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.trimmed.fq.gz ./HAM_1_EKDL240001890-1A_222TKYLT4_L1_2.orhpan.fq.gz \
        SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

```

# Stage 2: Mapping trimmed reads onto the C_excelsa_V5.fa reference genome with BWA

The Burrow-Wheeler aligner (BWA - version 12.3) was used to align trimmed reads onto the C_excelsa_V5.fa reference genome. An example command for BWA alignment can be found below.
```
#make environmental variables to store the metadata (meta), output directory (OUTDIR), and the path to the reference genome (REFDIR).
meta=EKDL240001890-1A_222TKYLT4
OUTDIR=~/220524_alignments
REFDIR=~/C_excelsa_V5_reference/C_excelsa_V5.fa

##Align the PAR_2 reads to C_excelsa_V5.fa reference genome
bwa mem \
     -t 16 \
      $REFDIR \
     ./PAR_2_${meta}_L1_1.trimmed.fq.gz ./PAR_2_${metad}_L1_2.trimmed.fq.gz \
     > $OUTDIR/PAR_2_${meta}_aln-pe.sam

```

# Stage 3: Converting to bam files, sorting, and indexing bams with Samtools

Samtools (version - 1.8) was used to convert the `aln-pe.sam` files produced by BWA in the previous stage into binary bam files using the `samtools view` command. The bam files were coordinate sorted using `samtools sort`, and subsequently indexed with `samtools index`. Finally, summary statistics for the alignments were produced with `samtools flagstats`. 

```
#samtools converting PAR_2 sam alignment to bam with samtools view 
samtools view -@ 4 -h -b ./PAR_2_${metadata}_aln-pe.sam -o ./bam_files/PAR_2_${metadata}.bam

#samtools sort to produce a coordinate-sorted bam file
samtools sort -@ 4 -o ./bam_files/PAR_2_${metadata}.sorted.bam ./bam_files/PAR_2_${metadata}.bam

#samtools index to produce an indexed sorted.bam file
samtools index ./bam_files/PAR_2_${metadata}.sorted.bam

#samtools flagstat produces alignment summary statistics
samtools flagstat ./bam_files/PAR_2_${metadata}.sorted.bam > ./bam_files/PAR_2_${metadata}.flagstats

```

# Stage 4: Marking and discarding duplicate reads with Picard MarkDuplicates and Adding Read Groups with AddOrReplaceReadGroups

Duplicate reads from the sorted bams were marked and discarded using Picard (version 3.0.0) MarkDuplicates, specifying `--REMOVE_DUPLICATES true`. An example MarkDuplicates command for one of the samples can be found below.

```
##make environmental variables for the output directory (OUTDIR) and the metadata (meta)
OUTDIR=~/220524_alignments/bam_files/duplicate_marked_bams
meta=EKDL240001890-1A_222TKYLT4

##execute MarkDuplicates on FLEET_2
java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ./FLEET_2_${meta}.sorted.bam -O $OUTDIR/FLEET_2_${meta}.marked_duplicates.bam -M $OUTDIR/FLEET_2_${meta}.marked_dup_metrics.txt --VALIDATION_STRINGENCY SILENT --ASSUME_SORTED true --REMOVE_DUPLICATES true

```

Read Groups were manually added to the duplicate marked bam files and were coordinate-sorted and indexed using Picard AddOrReplaceReadGroups with the `SORT_ORDER=coordinate` and `CREATE_INDEX=True` command line options. An example command for adding read groups can be found below.
```
meta=EKDL240001890-1A_222TKYLT4
#################
##fix read groups
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
```

# Stage 5: Genotyping Individual Samples with GATK HaplotypeCaller

Samples were genotyped utilising GATK (version 4.4.0) HaplotypeCaller specifying `--emit-ref-confidence BP_RESOLUTION`, `--minimum-mapping-quality-score 25` and `--min-base-quality-score 25`. An example command for a single accession/population can be found below.

```
##make environmental variables for the reference genome and for the output directory
OUTDIR=~/300524_HaplotypeCaller_output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
###########
##HaplotypeCaller for FLEET_2
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./FLEET_2_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
   -O $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -bamout $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.6x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 6 \
```

# Stage 6: Combining per-sample gVCFs into a single gVCF with GATK CombineGVCFs

The per-sample gVCFs generated in the previous stage by GATK HaplotypeCaller were combined into a multi-sample gVCF with `GATK CombineGVCFs`. The command used to combine all the gVCFs into a multi-sample gVCF can be seen below.

```
###make environmental variables for the reference genome and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
OUTDIR=~/300524_HaplotypeCaller_output/090624_Combined_VCF
################
##GATK CombineGVCFs of the additional danica and ionopsidium samples
 gatk CombineGVCFs \
   -R $REF \
   --variant ./Iac.g.vcf.gz \
   --variant ./Ime.g.vcf.gz \
   --variant ./Pen_1_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   --variant ./NOT_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   --variant ./SPEY_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   --variant ./LWS_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   --variant ./Iab_1.g.vcf.gz \
   --variant ./Iab_2.g.vcf.gz \
   --variant ./FLEET_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   --variant ./PAR_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -O $OUTDIR/090624_Ionops_danica_combined.g.vcf.gz
################
```

# Stage 7: Joint genotyping with GATK GenotypeGVCFs

The multi-sample gVCF produced in the previous stage with GATK CombineGVCFs was joint-genotyped using GATK GenotypeGVCFs specifying `-G StandardAnnotation` and `--include-non-variant-sites True`. The command used to joint-genotype the multi-sample gVCF can be found below.

```
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
INDIR=~/300524_HaplotypeCaller_output/090624_Combined_VCF
OUTDIR=~/300524_HaplotypeCaller_output/090624_combined_genotyped
################
##GATK GenotypeGVCFs of the additional danica samples and ionopsidium samples
 gatk GenotypeGVCFs \
   -R $REF \
   -V $INDIR/090624_Ionops_danica_combined.g.vcf.gz \
   -O $OUTDIR/090624_Ionops_danica_genotypd.g.vcf.gz \
   -G StandardAnnotation \
   --include-non-variant-sites True
################
```

# Stage 8: Filtering with GATK SelectVariants and GATK VariantFiltration

GATK SelectVariants was used to exclude insertion-deletion mutations and mixed SNP-indels, and to include only biallelic SNPs from the multi-sample gVCF.

```
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/300524_HaplotypeCaller_output/090624_combined_genotyped/090624_Ionops_danica_genotyped.g.vcf.gz
OUT1=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F1.biallelic.g.vcf.gz
OUT1=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F2.best.practice.g.vcf.gz
################
##GATK SelectVariants to select biallelic variants only
gatk SelectVariants \
   -R $REF \
   -V $VCF \
   -O $OUT1 \
   --select-type-to-exclude INDEL \
   --select-type-to-exclude MIXED \
   --restrict-alleles-to BIALLELIC \
```

Finally, VCFtools (version 1.16) was used to output per-site depth statistics.

```
#####
##use vcftools to output per site depth statistics
vcftools --gzvcf $OUT1 --out Depth.per.site --site-depth
```

## SplitsTree

SplitsTree was downloaded following the instructions on the [University of Tübingen Website](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/). 

This program was used to construct and visualize phylogenetic networks of the individuals and the populations in the ld pruned VCF file. (UPDATED ON THE 13th MAY 2024 - ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz).

Editing the phylogenetic networks was performed using Microsoft Powerpoint and manually highlighting clades. 

The 140524_adegenet_VCFs.R script was used to analyse the LD pruned and filtered VCF, utilising the glPcaFast() and vcf2genlightTetra() functions provided by Yant et al (2023). The VCF is loaded into Rstudio and is subsequently converted into a genlight object using the vcf2genlightTetra() function for polyploid data. Next, principal component analysis (PCA) can be performed on the genlight object, and subsequently, the genlight object can be converted into Nei's genetic distances using the stamppNeisD() function. Nei's genetic distances can be calculated for both the individual samples and the populations, and can be subsequently prepared for exporting into SplitsTree by the stamppPhylip() function provided by former student Anna (INSERT SURNAME, YEAR).   

## IQTREE and iTOL for maximum likelihood tree estimation and visualization

[IQTREE](http://www.iqtree.org/#download) was downloaded locally following the download instructions for `64-bit macOS Universal`. 

After navigating to the directory where the IQTREE executable is located, the following command was executed.

```
##Execute iqtree2 using Nei's genetic distance data and 4 threads/CPUs
bin/iqtree2 -s ~/Desktop/110624_aa.indiv_Neis_distance_4ds.phy -nt 4
```
The .iqtree file produced as output suggested that the substitution model that produces the maximum likelihood tree was `MK+I{0.0727447}+G4{0.244664}`, therefore, the analysis was re-done, this time with 1000 Bootstrap replicates for estimating branch supports and utilising the `-bnni` flag to reduce the risk of over-estimating branch supports. The command can be found below.
```
##execute iqtree2 with 4 threads, 1000 bootstrap replicates, the ML substitution model, and -bnni to reduce the risk of branch-support overestimation
~/Desktop/iqtree-2.3.4-macOS/bin/iqtree2 -s ~/Desktop/110624_IQTREE.OUT/110624_aa.indiv_Neis_distance_4ds.phy -nt 4 -B 1000 -m "MK+I{0.0727447}+G4{0.244664}" -bnni -redo
```

[iTOL](https://itol.embl.de/upload.cgi) or the Interactive Tree of Life, is a GUI which was used to upload the Newick-formatted consensus tree produced by IQTREE and to visualize the consensus tree. 

## Twisst (Topology weighting by iterative sampling of sub-trees)

This software can be used to quantify relationships between taxa that are not necessarily monophyletic, and can be used to explore how relationships between taxa varies across the genome by using genomic single nucleotide polymorphism (SNP) windows.

Citation: Simon H Martin, Steven M Van Belleghem, Exploring Evolutionary Relationships Across the Genome Using Topology Weighting, Genetics, Volume 206, Issue 1, 1 May 2017, Pages 429–438, [https://doi.org/10.1534/genetics.116.194720](https://doi.org/10.1534/genetics.116.194720).

The twisst.py script was downloaded from [Simon Martin's github page](https://github.com/simonhmartin/twisst/blob/master/twisst.py) and the [ete3 toolkit](http://etetoolkit.org/download/) and [numpy](https://numpy.org/) were downloaded following instructions on their respective websites.

The following commands were used to download and install both ete3 and numpy locally:
```
##create the ete3 environment
conda create -n ete3 python=3
##activate the ete3 environment
conda activate ete3
##use pip to install ete3 since conda didn't work initially
pip install ete3
##also use pip to install numpy
pip install numpy
##because of the errors when executing ete3 build check, install the packages required
pip install six
pip install pyqt5
pip install lxml
##for the final installation 
conda install -c etetoolkit ete_toolchain
##to check the ete3 toolkit has been installed correctly
ete3 build check
```






