# Genomic-Swamping-and-Introgression
This repository should allow the user to reproduce an analysis of the extent of introgression between the invasive hexaploid Cochlearia danica and native species Cochlearia officianalis, utilising tools such as Dsuite and Twisst. 

# Background

The *Cochlearia* species complex is a promising study system for the evolutionary genomics of adaptation due to the rapid acquisition of traits such as cold tolerance, salt tolerance, and heavy-metal tolerance, in a short evolutionary timescale ([Wolf et al., 2021](https://doi.org/10.7554/eLife.71572). This species complex thrive in cold Alpine and Arctic environments which contrasts greatly from the preferred and only Mediterranean habitat of its sole sister taxa *Ionopsidium* (Wolf et al., 2021). The genus *Cochlearia* consists of 16 accepted species and 4 subspecies, and started to rapidly diversify due to periodic climatic fluctuations that occurred during the Middle and Late Pleistocene (between 0.77-0.012 million years ago; Wolf et al., 2021). The genus *Cochlearia* displays a wide range of cytotypes and levels of ploidy, ranging from diploids to tetraploids to hexaploids to octaploids, which makes this species complex an interesting model species for the study of adaptations to whole genome duplication ([Bray et al., 2020](https://www.biorxiv.org/content/10.1101/2020.03.31.017939v1.full)). 

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


