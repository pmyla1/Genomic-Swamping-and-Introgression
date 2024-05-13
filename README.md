# Genomic-Swamping-and-Introgression
This repository should allow the user to reproduce an analysis of the extent of introgression between the invasive hexaploid *Cochlearia danica* and native species *Cochlearia officianalis*, utilising tools such as Dsuite and Twisst. 

# Background

The *Cochlearia* species complex is a promising study system for the evolutionary genomics of adaptation due to the rapid acquisition of traits such as cold tolerance, salt tolerance, and heavy-metal tolerance, in a short evolutionary timescale ([Wolf et al., 2021](https://doi.org/10.7554/eLife.71572). This species complex thrive in cold Alpine and Arctic environments which contrasts greatly from the preferred and only Mediterranean habitat of its sole sister taxa *Ionopsidium* (Wolf et al., 2021). The genus *Cochlearia* consists of 16 accepted species and 4 subspecies, and started to rapidly diversify due to periodic climatic fluctuations that occurred during the Middle and Late Pleistocene (between 0.77-0.012 million years ago; Wolf et al., 2021).

The genus *Cochlearia* displays a wide range of cytotypes and levels of ploidy, ranging from diploids to tetraploids to hexaploids to octaploids, which makes this species complex an interesting model species for the study of adaptations to whole genome duplication ([Bray et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.017939v1.full)). Since whole genome duplication and polyploidy are major effect mutations, many cellular, ionomic, and molecular processes are disrupted, especially those pertaining to sister chromatid segregation during meiosis, DNA repair, and recombination, amongst others ([Yant & Schmickl, 2021](https://pubmed.ncbi.nlm.nih.gov/33454987/)). Therefore, neopolyploids such as some *Cochlearia* species must overcome these initial challenges associated with sister chromatid segregation during meiosis, to enable reduced crossover numbers and to prevent chromosomal breakages during anaphase (Bray et al., 2020). Some of the genes found to be under selection in neopolyploids relative to their diploid counterparts are involved in biological processes such as DNA repair, recombination, sister chromatid segregation, amongst others, however, despite this process-level convergence there appears to be low orthologue-level convergence ([Bray et al., 2023](https://www.biorxiv.org/content/10.1101/2023.09.27.559727v1.full)). 

*Cochlearia danica* is a highly invasive, salt tolerant species in the *Cochlearia* genus that is native to the Atlantic coasts of Europe in countries including Denmark, Belgium, the UK, Germany, Finland, Ireland, Norway, Russia, and Sweden amongst others (Fekete et al., 2018). Despite being native to the Atlantic coasts of Europe, *C. danica* has spread rapidly throughout Central Europe since the 1970s, especially along roadsides which is attributed to the widespread use of de-icing salts ([Fekete et al., 2018](http://dx.doi.org/10.23855/preslia.2018.023)). The rate of spread of *C. danica* along Central European roadsides was estimated to be approximately 62-65km/year, and the soil in which it grows and thrives is characterised by high salt content (Fekete et al., 2018). Interestingly, this species has been reported to undergo rapid and marked/remarkable changes in population size, and one Hungarian population was found to decrease in size by 99% between 2016 and 2017 (Fekete et al., 2018).


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

## Twisst (Topology weighting by iterative sampling of sub-trees)

This software can be used to quantify relationships between taxa that are not necessarily monophyletic, and can be used to explore how relationships between taxa varies across the genome by using genomic single nucleotide polymorphism (SNP) windows.

Citation: Simon H Martin, Steven M Van Belleghem, Exploring Evolutionary Relationships Across the Genome Using Topology Weighting, Genetics, Volume 206, Issue 1, 1 May 2017, Pages 429–438, [https://doi.org/10.1534/genetics.116.194720](https://doi.org/10.1534/genetics.116.194720).


## SplitsTree

SplitsTree was downloaded following the instructions on the [University of Tübingen Website](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/). 

This program was used to construct and visualize phylogenetic networks of the individuals and the populations in the ld pruned VCF file.

Editing the phylogenetic networks was performed using Microsoft Powerpoint and manually highlighting clades.  

## Genome Analysis Toolkit (GATK)

## R & Rstudio








