# Genomic-Swamping-and-Introgression
This repository should allow the user to reproduce an analysis of the extent of **introgression** between the invasive hexaploid ***Cochlearia danica*** and native tetraploid species ***Cochlearia officianalis***, utilising tools such as **Dsuite** and **Twisst**. 

# Background

The ***Cochlearia* species complex** is a promising study system for the **evolutionary genomics of adaptation** due to the rapid acquisition of traits such as **cold tolerance**, **salt tolerance**, and **heavy-metal tolerance**, in a short evolutionary timescale ([Wolf et al., 2021](https://doi.org/10.7554/eLife.71572)). This species complex **thrives** in cold **Alpine and Arctic** environments which contrasts greatly from the preferred **Mediterranean habitat** of its sole **sister taxa *Ionopsidium*** (Wolf et al., 2021). The genus *Cochlearia* consists of **16 accepted species** and **4 subspecies**, and started to **rapidly diversify** due to periodic **climatic fluctuations** that occurred during the **Middle** and **Late Pleistocene** (between **0.77-0.012** million years ago; Wolf et al., 2021).

The genus *Cochlearia* displays a **wide range of cytotypes** and levels of ploidy, ranging from **diploids** to tetraploids to hexaploids to **octaploids**, which makes this species complex an interesting **model species** for the study of **adaptations to whole genome duplication** ([Bray et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.017939v1.full)). Since whole genome duplication and polyploidy are **major effect mutations**, many **cellular**, **ionomic**, and **molecular** processes are **disrupted**, especially those pertaining to **sister chromatid segregation** during meiosis, **DNA repair**, and **recombination**, amongst others ([Yant & Schmickl, 2021](https://pubmed.ncbi.nlm.nih.gov/33454987/)). Therefore, **neopolyploids** must overcome these initial challenges associated with sister chromatid segregation during meiosis, to enable **reduced crossover numbers** and to **prevent chromosomal breakages** during anaphase (Bray et al., 2020). Some of the **genes** found to be **under selection in neopolyploids** relative to their diploid counterparts are involved in biological processes such as **DNA repair**, **recombination**, **sister chromatid segregation**, amongst others, however, despite this **process-level convergence** there appears to be **low orthologue-level convergenc** ([Bray et al., 2023](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v1.full)). 

*Cochlearia danica* is a **highly invasive, salt tolerant** species in the *Cochlearia* genus that is **native** to the **Atlantic coasts of Europe** in countries including Denmark, Belgium, the UK, Germany, Finland, Ireland, Norway, Russia, and Sweden amongst others (Fekete et al., 2018). Despite being native to the Atlantic coasts of Europe, *C. danica* has **spread rapidly throughout Central Europe** since the 1970s, especially **along roadsides** which is attributed to the **widespread use of de-icing salts** ([Fekete et al., 2018](http://dx.doi.org/10.23855/preslia.2018.023)). The **rate of spread** of *C. danica* along Central European roadsides was estimated to be **approximately 62-65km/year**, and the soil in which it grows and thrives is characterised by **high salt content** (Fekete et al., 2018). Interestingly, this species has been reported to undergo **rapid and marked/remarkable** changes in **population size**, and one Hungarian population was found to decrease in size by 99% between 2016 and 2017 (Fekete et al., 2018).


# Installation of Software and Dependencies

## Dsuite
[Dsuite](https://github.com/millanek/Dsuite) is a software program developed to quickly calculate Patterson's D (ABBA-BABA), and the f4-ratio statistics across many populations and/or species. 

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

The TREE_FILE.nwk was structured as following, where the Outgroup was the UK diploid species *Cochlearia pyrenaica*:
```
(Outgroup,(C_officinalis,(C_anglica,C_danica)));
```
The SETs.txt file was structured following the guidelines on the [Dsuite](https://github.com/millanek/Dsuite) Github page with the individual ID and the group ID separated by a tab, and is demonstrated below:
```
BNK21  Outgroup
CHA_1  Outgroup
CHA_2  Outgroup
JOR_1  Outgroup
...
AAH_1  C_officinalis
AAH_2  C_officinalis
AAH_3  C_officinalis
AAH_4  C_officinalis
...
BRE_1  C_danica
CUM_1  C_danica
DAR_1  C_danica
DAR_3  C_danica
...
SKF_002  C_anglica
SKF_003  C_anglica
SKF_005  C_anglica
SKF_009  C_anglica
```

## 15/05/2024 Additional Dsuite commands and calculations

Dsuite Dtrios was also executed on a newly LD pruned, filtered VCF file containing all UK hexaploids, including all *C. danica* and putative *C. anglica* samples, all UK diploids, some EU diploids, and approximately half of the UK and EU tetraploids.

This time, the TREE_FILE.nwk was altered to change the relationships between the sequences

```
##utilise Dsuite
cd ~/Dsuite/

##execute Dtrios with a jackknife number of 50 (splits the SNPs into 50 blocks)
./Build/Dsuite Dtrios -k 50 -t ../Desktop/Individual_Project_files/Dsuite_files/150524_experimental_TREE_FILE.nwk --ABBAclustering -o 150524_EU_UK_ld_pruned ../Desktop/Individual_Project_files/VCF_files/150524_ld_pruned_20PCTmis_maf005_EU_UK_pops.vcf.gz ../Desktop/Individual_project_files/Dsuite_files/150524_SETs.txt
```

The structure of the 150524_experimental_TREE_FILE.nwk can be found below, where *C. danica* replaces *C. officinalis*:

```
(Outgroup,(C_danica,(C_anglica,C_officinalis)));
```

# 17/05/24 

To extract the individual samples names from the VCF file use bcftools and the following command:
```
bcftools query -l 160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex.vcf.gz > ./170524_allUKtets_allUKhets_samples.txt
```
# Fastqc and Multiqc for sequencing quality control reports

The Multiqc report for the additional *Cochlearia danica* and *Ionopsidium* sequencing data provided by Yant (2024) can be accessed via the following link [INSERT LINK HERE]. 

An example command for producing a fastqc sequencing quality control report can be found below:
```
###load the fastqc module for sequencing quality control
module load fastqc-uoneasy/0.12.1-Java-11

##change directory to where your fastq.gz files are stored
cd /file/path/to/fastq.gz files

##perform fastqc on the FLE_2 sequencing data/reads
fastqc -o ../170524_fastqc/ ./FLE_2/*.fq.gz

##repeat for the other populations

##unload fastqc module
module unload fastqc-uoneasy/0.12.1-Java-11
```

Subsequently, a multiqc report can be performed/executed on the directory containing the results from the fastqc reports, which will be .fastqc.gz files. The command required to produce the multiqc report can be found below.
```
##load multiqc module
module load multiqc-uoneasy/1.14-foss-2023a

##execute multiqc on the fastqc.zip data specifying -f (--force to overwrite existing reports) and -p to export the plots generated
multiqc /gpfs01/home/pmyla1/170524_fastqc/.*fastqc.zip

##unload multiqc module
module unload multiqc-uoneasy/1.14-foss-2023a
```
The multiqc plots and reports can be found within the directory from which the commands were executed and include a directory for png, svg, or pdf versions of the plots.

# Cutadapt (Trimming adapters from Illumina paired-end reads)

Cutadapt and the 170524_cutadapters_fastqc_multiqc.sh script was used to trim the Nextera transposase adapter sequences from the fastq.gz sequencing files. 

Example command:
```
##cut adapters from the end of the fastq.gz files with cutadapt
##the read 1 Nextera transposase sequence is TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
## the read 2 Nextera transposase adapter sequence is GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
module load cutadapt-uon/gcc12.3.0/4.6

##first cut the nextera transposase read 1 sequence from FLE_2
cutadapt -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz ./FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz

##then cut the nextera transposase read 2 sequence from FLE_2
cutadapt -a GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o 170524_cutadapt/170524_FLEET_2_EKDL240001890-1A_222TKYLT4_L1_2.fq.gz ./FLE_2/FLEET_2_EKDL240001890-1A_222TKYLT4_L1_1.fq.gz

```

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

## SplitsTree

SplitsTree was downloaded following the instructions on the [University of Tübingen Website](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/). 

This program was used to construct and visualize phylogenetic networks of the individuals and the populations in the ld pruned VCF file. (UPDATED ON THE 13th MAY 2024 - ld_pruned_130524_allUKdips_allUKtets_allUKhex.vcf.gz).

Editing the phylogenetic networks was performed using Microsoft Powerpoint and manually highlighting clades. 

The 140524_adegenet_VCFs.R script was used to analyse the LD pruned and filtered VCF, utilising the glPcaFast() and vcf2genlightTetra() functions provided by Yant et al (2023). The VCF is loaded into Rstudio and is subsequently converted into a genlight object using the vcf2genlightTetra() function for polyploid data. Next, principal component analysis (PCA) can be performed on the genlight object, and subsequently, the genlight object can be converted into Nei's genetic distances using the stamppNeisD() function. Nei's genetic distances can be calculated for both the individual samples and the populations, and can be subsequently prepared for exporting into SplitsTree by the stamppPhylip() function provided by former student Anna (INSERT SURNAME, YEAR).   

## RAXML 

Raxml was downloaded following the instructions on the [raxml-ng](https://github.com/amkozlov/raxml-ng?tab=readme-ov-file) Github page, by selecting the **Download OSX/macOS binary (x86 and ARM)** link on the page.

To display the help messages you can execute the following command from the directory of the raxml-ng executable:
```
./raxml-ng -h 
```

## Genome Analysis Toolkit (GATK)

## R & Rstudio








