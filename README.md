# Genomic-Swamping-and-Introgression
This repository should allow the user to reproduce an analysis of **gene flow/introgression** between the rapidly spreading invasive hexaploid ***Cochlearia danica*** and native UK ***Cochlearia species*** including: the diploid ***Cochlearia pyrenaica***, the tetraploid ***Cochlearia officinalis***, and the hexaploid ***Cochlearia anglica***-like hybrid, utilising the **Dsuite** software package for calculating Patterson's D and f4-admixture ratios (ABBA-BABA statistics). 

# Background

The ***Cochlearia* species complex** is a promising study system for the **evolutionary genomics of adaptation** due to the rapid acquisition of traits such as **cold tolerance**, **salt tolerance**, and **heavy-metal tolerance**, in a short evolutionary timescale ([Wolf et al., 2021](https://doi.org/10.7554/eLife.71572)). This species complex **thrives** in cold **Alpine and Arctic** environments in contrast to the preferred **Mediterranean habitat** of its **sister taxa *Ionopsidium*** (Wolf et al., 2021). The genus *Cochlearia* consists of **16 accepted species** and **4 subspecies**, and started to **rapidly diversify** due to periodic **climatic fluctuations** during the **Middle** and **Late Pleistocene** (between **0.77-0.012** million years ago; Wolf et al., 2021).

The genus *Cochlearia* displays a **wide range of cytotypes** and ploidies, ranging from **diploids** through **autotetraploids** to **hexaploids** and **octaploids**, which makes *Cochlearia* an interesting **model system** for studying **adaptation to whole genome duplication (WGD)** ([Bray et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.017939v1.full)). WGD and polyploidy are **major effect mutations** which drastically disrupt **cellular**, **ionomic**, and **molecular** processes, especially relating to **sister chromatid segregation** during meiosis, **DNA repair**, and **recombination** ([Yant & Schmickl, 2021](https://pubmed.ncbi.nlm.nih.gov/33454987/)). Therefore, newly formed polyploids or **neo-polyploids** must overcome the initial challenges associated with sister chromatid segregation during meiosis, ultimately to **prevent chromosomal breakages** during anaphase (Bray et al., 2020). Some of **candidate genes under selection in neo-polyploids** function in biological processes like **DNA repair**, **recombination**, **sister chromatid segregation**, however, despite **process-level convergence** there is **low orthologue-level convergence** ([Bray et al., 2023](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v1.full)). 

*Cochlearia danica* is a **highly invasive, salt tolerant** species in the *Cochlearia* genus that is **native** to the **Atlantic coasts of Europe** but has rapidly spread throughout **continental Europe since the 1970s** which has been attributed to the **widespread use of de-icing salts** ([Fekete et al., 2018](http://dx.doi.org/10.23855/preslia.2018.023)). The **rate of spread** of *C. danica* along Central European roadsides was estimated to be **approximately 62-65km/year**, *C. danica* may experience **rapid and marked** changes in **population size**, which is emphasised by the **99% reduction in population size** of one Hungarian population between **2016 and 2017** (Fekete et al., 2018). The rapid population dynamics of *C. danica* have led to the hypothesis that this invasive halophyte could be breeding with native UK octaploid *Cochlearia anglica* resulting in **genomic swamping**. Supporting this hypothesis is the phenotypic observation of hexaploid *Cochlearia anglica*-like hybrids in the Eastern United Kingdom which suggests that there is introgression between *C. danica* and native UK accessions of the octaploid *C. anglica*. 

**Introgression** is the exchange of genetic material between species that results from **hybridization** and **recurrent backcrossing** ([Wang et al., 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10504873/)). Introgression has the potential to introduce large sets of new alleles simultaneously at multiple unlinked loci, which permits adaptation even in polygenic traits thus potentially promoting rapid species evolution ([Burgarella et al., 2019](https://doi.org/10.3389/fpls.2019.00004)). Introgression is considered **"adaptive"** if the genetic material transferred confers the **recipient species an increased fitness**, which occurs between **crop wild relatives** and their **domesticated counterparts** providing the latter with beneficial traits such as **increased resistance** to **biotic and abiotic stress** including drought or extreme temperatures ([Burgarella et al., 2019](https://doi.org/10.3389/fpls.2019.00004)).

The aim of this project is to detect and quantify the extent of gene flow/introgression between the invasive hexaploid halophyte *C. danica* and native UK populations of *C. anglica* (6X) or UK populations of *C. officinalis* (4X) by calculating ABBA-BABA statistics on genome-wide SNP data using Dsuite ([Malinsky, 2021](https://doi.org/10.1111/1755-0998.13265)). 


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

# Data Generation

Additional *C. danica* and *Ionopsidium* Illumina paired-end sequencing data was provided by Yant (2024), and the sequencing reads (in fq.gz/fastq.gz format) were processed following the guidelines outlined in [ngs_pipe](https://github.com/mattheatley/ngs_pipe) from Healey (2024).

## FastQC and MultiQC for sequencing quality control reports

An example command for producing a fastqc report for a single population (FLEET_2 in this example) can be found below:
```
###load fastqc module for performing sequencing quality control reports
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
multiqc -f -p ~/170524_fastqc/.*fastqc.zip
```
The MultiQC plots and reports include a directory for **png**, **svg**, or **pdf** versions of the plots.

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

The Burrows-Wheeler aligner (BWA - version 12.3) was used to align trimmed reads onto the C_excelsa_V5.fa reference genome. An example command for BWA alignment can be found below.
```
#make environmental variables to store the metadata (meta), output directory (OUTDIR), and the path to the reference genome (REFDIR).
meta=EKDL240001890-1A_222TKYLT4
OUTDIR=~/220524_alignments
REFDIR=~/C_excelsa_V5_reference/C_excelsa_V5.fa

##Align the PAR_2 reads to C_excelsa_V5.fa reference genome using 16 threads 
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
##store the metadata in an environmental variable called meta
meta=EKDL240001890-1A_222TKYLT4
##fix read groups with picard AddOrReplaceReadGroups
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
```

# Stage 5: Genotyping Individual Samples with GATK HaplotypeCaller

Samples were genotyped utilising GATK (version 4.4.0) HaplotypeCaller specifying `--emit-ref-confidence BP_RESOLUTION`, `--minimum-mapping-quality-score 25`, `--min-base-quality-score 25`, and `--sample-ploidy 6` because the additional *C. danica* samples are hexaploid. An example command for a **single population** can be found below.

```
##make environmental variables for the reference genome and for the output directory
OUTDIR=~/300524_HaplotypeCaller_output
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa

##HaplotypeCaller for FLEET_2
 gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R $REF \
   -I ./FLEET_2_EKDL240001890-1A_222TKYLT4.marked_duplicates.bam \
   -O $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.g.vcf.gz \
   -bamout $OUTDIR/FLEET_2_EKDL240001890-1A_222TKYLT4.6x.bam \
   --emit-ref-confidence BP_RESOLUTION \
   --min-base-quality-score 25 \
   --minimum-mapping-quality 25 \
   --sample-ploidy 6 
```

# Stage 6: Combining per-sample gVCFs into a single gVCF with GATK CombineGVCFs

The per-sample gVCFs generated in the previous stage by GATK HaplotypeCaller were combined into a **multi-sample gVCF** with `GATK CombineGVCFs`. The command used to combine all **per-sample gVCFs** into a **multi-sample gVCF** is below.

```
###make environmental variables for the reference genome and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
OUTDIR=~/300524_HaplotypeCaller_output/090624_Combined_VCF

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
```

# Stage 7: Joint genotyping with GATK GenotypeGVCFs

The multi-sample gVCF produced in the previous stage with GATK CombineGVCFs was **joint-genotyped** using GATK GenotypeGVCFs specifying `-G StandardAnnotation` and `--include-non-variant-sites True`. The command used to joint-genotype the multi-sample gVCF can be found below.

```
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
INDIR=~/300524_HaplotypeCaller_output/090624_Combined_VCF
OUTDIR=~/300524_HaplotypeCaller_output/090624_combined_genotyped

##GATK GenotypeGVCFs of the additional danica samples and ionopsidium samples
 gatk GenotypeGVCFs \
   -R $REF \
   -V $INDIR/090624_Ionops_danica_combined.g.vcf.gz \
   -O $OUTDIR/090624_Ionops_danica_genotypd.g.vcf.gz \
   -G StandardAnnotation \
   --include-non-variant-sites True
```

# Stage 8: Filtering with GATK SelectVariants and GATK VariantFiltration

GATK SelectVariants was used to **exclude insertion-deletion** mutations and **mixed SNP-indels**, and to **include only biallelic SNPs** from the multi-sample gVCF.

```
###make environmental variables for the reference genome, input directory, and the output directory
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/300524_HaplotypeCaller_output/090624_combined_genotyped/090624_Ionops_danica_genotyped.g.vcf.gz
OUT1=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F1.biallelic.g.vcf.gz
OUT1=~/300524_HaplotypeCaller_output/090624_combined_genotyped/100624_filtered.best/100624_Ion.dan.filtered.F2.best.practice.g.vcf.gz

##GATK SelectVariants to select biallelic variants only
gatk SelectVariants \
   -R $REF \
   -V $VCF \
   -O $OUT1 \
   --select-type-to-exclude INDEL \
   --select-type-to-exclude MIXED \
   --restrict-alleles-to BIALLELIC \
```

Finally, VCFtools (version 1.16) was used to output **per-site depth statistics**.

```
##use vcftools to output per site depth statistics and use 1.6 * mean depth as a maximum depth cut-off
vcftools --gzvcf $OUT1 --out Depth.per.site --site-depth
```

# Stage 9: Depth filtering and final VCF generation

GATK (version 4.4.0) SelectVariants was used to produce a depth-masked VCF file (120624_depth.mask.Ion.dan.g.vcf.gz.) based on a depth cut off of 1.6 * the mean depth. Subsequently, GATK VariantFiltration was used to filter the F2 best practice VCF (110624_Ion.dan.F2.best.g.vcf.gz) using the depth-masked VCF and removing the variants NOT in the masked VCF. The final F4 VCF (120624_Ion.dan.filtered.F4.g.vcf.gz) was produced using GATK SelectVariants by excluding the depth-filtered sites from the F3 VCF with the `--exclude-filtered True` command line option.
 

```
##environmental variables for reference genome (REF), input VCF (VCF), and output depth mask VCF (OUTMASK)
REF=~/C_excelsa_V5_reference/C_excelsa_V5.fa
VCF=~/300524_HaplotypeCaller_output/110624_Ion.dan.F2.best.g.vcf.gz
OUTMASK=~/300524_HaplotypeCaller_output/120624_depth.mask.Ion.dan.g.vcf.gz

##Use GATK SelectVariants to filter based on a maximum depth cut off of 1.6 * mean depth
gatk SelectVariants \
        -R $REF \
        -V $VCF \
        -O $OUTMASK \ ##produces depth-mask VCF 
        --select "DP<149" ##depth cut-off = 1.6 * mean depth

##environmental variable for the F3-depth filtered VCF
OUTF3=~/300524_HaplotypeCaller_output/120624_Ion.dan.filtered.F3.g.vcf.gz

gatk VariantFiltration \
        -R $REF \
        -V $VCF \
        -O $OUTF3 \
        --mask $OUTMASK \
        --filter-not-in-mask ##remove the sites NOT in the masked VCF

##environmental variable for final F4-filtered VCF
OUTF4=~/300524_HaplotypeCaller_output/120624_Ion.dan.filtered.F4.g.vcf.gz

gatk SelectVariants \
    -R $REF \
    -V $OUTF3 \ ##input F3-filtered VCF from previous stage
    -O $OUTF4 \ ##final output F4-VCF
    --exclude-filtered True ##exclude the sites with depth > 149x
```

# Stage 10: Manually re-headering and merging VCF files 

BCFtools (version 1.18) was used to re-header the original 133 sample reheadered.F4_133.ann.vcf.gz VCF and the 120624_Ion.dan.filtered.F4.g.vcf.gz containing the additional *C. danica* and *Ionopsidium* samples. Firstly, the original VCF file headers were visualised using the following commands:

```
##environmental variable for reheadered.F4_133.ann.vcf.gz
VCF1ORIGINAL=reheadered.F4_133.ann.vcf.gz
bcftools view -h $VCF1ORIGINAL > ./HEADER1.txt

##environmental variable for 120624_Ion.dan.filtered.F4.g.vcf.gz
VCF2ORIGINAL=120624_Ion.dan.filtered.F4.g.vcf.gz
bcftools view -h $VCF2ORIGINAL > ./HEADER2.txt


```
Next, the HEADER1.txt and HEADER2.txt files were manually altered using `nano` (or an equivalent text editor) and changing the Number in the PL field in the VCF header to `Number=.`. Then, `bcftools reheader` was used to reheader the VCFs to allow for merging.

```
##environmental variable for the reheadered VCFs
VCF1REHEADERED=120624_reheadered.F4_133.ann.vcf.gz
bcftools reheader -h HEADER1.txt $VCF1ORIGINAL > ./$VCF1REHEADERED

VCF2REHEADERED=120624_reheadered.Ion.dan.F4.g.vcf.gz
bcftools reheader -h HEADER2.txt $VCF2ORIGINAL > ./$VCF2REHEADERED
```

The reheadered VCFs were indexed using a simple `bcftools index` command, and the 120624_reheadered.F4_133.ann.vcf.gz and 120624_reheadered.Ion.dan.F4.g.vcf.gz files were merged using  `bcftools merge`.

```
##make an environmental variable for the merged output VCF
FINALMERGEDVCF=120624.final.merged.Ion.dan.F4.vcf.gz
##use bcftools merge with 8 threads to combine/merge the reheadered VCFs together, assuming missing genotypes are unphased (0/0)
bcftools merge --t 8 -0 --write-index -Oz $VCF1REHEADERED $VCF2REHEADERED -o $FINALMERGEDVCF
```

# Stage 11: Selecting only UK diploids, tetraploids, and hexaploids and Ionopsidium samples

The 120624.final.merged.Ion.dan.F4.vcf.gz was indexed using `gatk IndexFeatureFile` and all biallelic SNPs from the UK samples were selected using `gatk SelectVariants` with the `--select-type-to-include SNP` and `--restrict-alleles-to BIALLELIC` flags. The `-sn` flag was used to select the individual IDs for the UK samples only (e.g. AAH_1).  

```
INVCF=~/120624_merged.133.ann.Ion.dan.F4.vcf.gz
OUTVCF=~/300524_HaplotypeCaller_output/120624_Ionops.allUKsamples.F4.vcf.gz

##index the 120624.final.merged.Ion.dan.F4.vcf.gz file with gatk IndexFeatureFile
gatk IndexFeatureFile -I $INVCF

##use GATK SelectVariants to select all Ionopsidium, UK diploids, tetraploids, C. danica, and putative C. anglica
gatk SelectVariants -V $INVCF \
 --select-type-to-include SNP \
 --restrict-alleles-to BIALLELIC \
 -sn AAH_1 -sn AAH_2 -sn AAH_3 -sn AAH_4 \
 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 \
 -sn BNK_21 -sn BRE_1 -sn CHA_1 -sn CHA_2 -sn CUM_1 \
 -sn DAR_1 -sn DAR_3 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 \
 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FOR_1 -sn FRE_013 \
 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn GEO_2 -sn GEO_6 \
 -sn Ime -sn Iac -sn Iab_1 -sn Iab_2 \
 -sn JON_001 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 \
 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 \
 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 \
 -sn LOS_1 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 \
 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 \
 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 \
 -sn PAR_2 -sn Pen_1 -sn NOT -sn LWS \
 -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn RYE_1 -sn SPEY_2 -sn SCO_1 \
 -sn SCU_1 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 \
 -sn SKF_002 -sn SKF_003 -sn SKF_005 -sn SKF_009 \
 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 \
 -sn SPU_006 -sn SPU_008 -sn SPU_009 -sn SPU_010 \
 -sn TET_002 -sn TET_004 -sn TET_006 -sn TET_008 \
 -O $OUTVCF \
 --allow-nonoverlapping-command-line-samples
```

# Stage 12: LD pruning and final VCF generation

The [prune_ld.c](https://github.com/thamala/polySV/blob/main/prune_ld.c) script from Hämälä (2024) was used to thin the VCF and remove all sites with more than 10% missing data, a minor allele frequency less than 0.05, and a squared genotypic correlation of 0.1 (in windows of 50 SNPs and a window step size of 10 SNPs). 

The script was first compiled using the following command:

```
##compile the prune_ld.c script 
gcc prune_ld.c -o prune_ld -lm 
```

The 120624_Ionops.allUKsamples.F4.vcf.gz was first unzipped using `gunzip` from the samtools suite and subsequently LD-pruned using the prune_ld executable.

```
##LD-prune the VCF - maximum missing 10%, minor allele frequency 0.05, squared genotypic correlation 50 10 0.1
prune_ld -vcf ~/120624_Ionops.allUKsamples.F4.vcf -mis 0.9 -maf 0.05 -r2 50 10 0.1 > ~/120624.LD.Pruned.Ionops.allUKsamples.vcf
```

# Phylogenetic Trees and Principal Component Analysis

## SplitsTree

SplitsTree was downloaded following the instructions on the [University of Tübingen Website](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/). 

This program was used to construct and visualize phylogenetic networks of the individuals in the ld pruned 120624.LD.Pruned.Ionops.allUKsamples.vcf. 

Editing the phylogenetic networks was performed using Microsoft Word and manually highlighting clades according to ploidy. 

The 030724.adegenet.R script was used to analyse the LD pruned and filtered VCF (120624.LD.Pruned.Ionops.allUKsamples.vcf), utilising the glPcaFast() and vcf2genlightTetra() functions provided by Yant et al (2023). 

The VCF was loaded into Rstudio and converted into a genlight object using the vcf2genlightTetra() function for polyploid data. Next, principal component analysis (PCA) can be performed on the genlight object using the glPcaFast() function, and subsequently, the genlight object can be converted into Nei's genetic distances using the stamppNeisD() function. 

Nei's genetic distances can be calculated for both the individual samples and the populations, and can be subsequently prepared for exporting into SplitsTree by the stamppPhylip() function.   

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
bin/iqtree2 -s ~/Desktop/110624_IQTREE.OUT/110624_aa.indiv_Neis_distance_4ds.phy -nt 4 -B 1000 -m "MK+I{0.0727447}+G4{0.244664}" -bnni -redo
```

[iTOL](https://itol.embl.de/upload.cgi) or the Interactive Tree of Life, is a GUI which was used to upload the Newick-formatted consensus tree produced by IQTREE and to visualize the consensus tree. 

To visualise your consensus tree you can upload the consensus tree in Newick format into the `Tree Text` box and select upload. Next you can customise the layout of your consensus tree as you wish by selecting the toolbar which includes `Basic`, `Advanced`, and `Datasets`. 

# Dsuite : Fast ABBA-BABA statistics and F4-admixture ratio calculations


## Dquartets - a programme to detect introgression between a quartet of species without an outgroup

Dquartets is part of the Dsuite software package from [Malinsky, 2021](https://github.com/millanek/Dsuite), and can be used to calculate the ABBA-BABA and F4-admixture ratio statistics for all possible quartets of species and does not require an outgroup. The species in the `SETS_SPECIES.txt` file were the individual IDs (3 letter population code followed by a number, e.g. AAH_1) and the species ID (*pyrenaica*, *officinalis*, *anglica*, or *danica*) separated by a tab. 

```
##SETS_SPECIES.txt file format
BNK_21       pyrenaica
...
AAH_1        officinalis
...
SKF_002      anglica
...
BRE_1        danica

```
Dquartets was executed on the 120624_LD.Pruned.Ionops.allUKsamples.vcf.gz using a jack-knife block approach to split the VCF into 4000 blocks of approximately 92 SNPs each (366,504 SNPs in total).  

```
##execute Dquartets on the dataset to obtain the assumed relationships between the species excluding the Ionopsidium outgroup samples
Dsuite Dquartets -k 4000 -o 120624_experimental $VCF SETS_SPECIES.txt
```

## Dsuite Dtrios - a fast programme to calculate ABBA-BABA and F4-admixture ratio statistics

Dtrios was used to calculate Patterson's D (ABBA-BABA) and F4-ratio statistics for all possible trios of species using *Ionopsidium* as an outgroup in the analysis. 

The `--ABBAclustering` option was used to test whether strong ABBA-informative sites cluster together throughout the genome. If introgression has occurred between two species, you would expect clusters of ABBA-informative sites across the genome rather than having many individual ABBA-informative sites evenly distributed across the genome caused by homoplasy, therefore, the `--ABBAclustering` option can be used to test for clustering of ABBA-informative sites. The more significant clustering of ABBA sites, the more confidence you can have that the introgression/gene flow event is real and not caused by homoplasies ([Malinsky, 2021](https://github.com/millanek/Dsuite)).

In order to execute Dsuite commands locally (e.g. Dtrios), you can navigate to the Build folder and run the Dsuite executable with the following command `./Build/Dsuite` which shows the available commands. To execute the Dtrios command you can type `./Build/Dsuite Dtrios`.

Dsuite Dtrios was executed using a modular approach on Ada. Dtrios was executed using a jack-knife block approach which divides the 120624_LD.Pruned.Ionops.allUKsamples.vcf.gz file into 4000 blocks of approximately 92 single nucleotide polymorphisms (total number of biallelic SNPs in the 120624_LD.Pruned.Ionops.allUKsamples.vcf.gz file is  366,504):

```
##make an environmental variable for the VCF you want to use (120624_LD.Pruned.Ionops.allUKsamples.vcf.gz)
VCF=~/120624_LD.Pruned.Ionops.allUKsamples.vcf.gz

#execute Dtrios using the appropriate outgroup (Iac, Ime, Iab_1, Iab_2) using 4000 Jack-knife blocks (-k 4000)
#without specifying an explictly stated phylogenetic tree
Dsuite Dtrios -k 4000 -o 240624_Dtrios --ABBAclustering $VCF SETS_SPECIES.txt
```

The SETs.txt file has the following structure with the individual ID and the group/species ID (i.e. the species) separated by a tab, and is demonstrated below:

```
Ime          Outgroup
Iac          Outgroup
Iab_1        Outgroup
Iab_2        Outgroup
...         
BNK_21       pyrenaica
CHA_1        pyrenaica
CHA_2        pyrenaica
JOR_1        pyrenaica
...
AAH_1        officinalis
AAH_2        officinalis
AAH_3        officinalis
AAH_4        officinalis
...
BRE_1        danica
CUM_1        danica
DAR_1        danica
DAR_3        danica
...
SKF_002      anglica
SKF_003      anglica
SKF_005      anglica
SKF_009      anglica
```

## Dinvestigate - a window-based introgression scan in trios with significantly elevated D-statistics

Dinvestigate was used to perform a window-based scan for introgression in trios that had significantly elevated D-statistics from the Dtrios output. 

The `240624_testtrios.txt` is a text file containing the trio of populations/species to test for localised regions of introgression separated by a tab:
```
##240624_testtrios.txt structure
pyrenaica        officinalis        anglica
```

An example command for one of the trios with elevated D-statistics (*C. pyrenaica*        *C. officinalis*        *C. anglica*) is shown below with two different SNP window sizes.
```
##first calculate windowed D-statistics for pyrenaica        officinalis        anglica trio (50 SNPs, 25 SNP step size)
Dsuite Dinvestigate -w 50,25 -n 50_25_pyr_off_ang $VCF SETS_SPECIES.txt 240624_testtrios.txt
##now calculate windowed D-statistics for pyrenaica        officinalis        anglica trio (100 SNPs, 25 SNP step size)
Dsuite Dinvestigate -w 100,25 -n 100_25_pyr_off_ang $VCF SETS_SPECIES.txt 240624_testtrios.txt
```

The output for the Dinvestigate analysis includes text files containing the localised windowed F-statistics including Fd, Fdm, and df. For subsequent/downstream analyses, the text files can be uploaded into RStudio and the top 1% introgression windows (using f_dM) can be selected using a dplyr based approach.



