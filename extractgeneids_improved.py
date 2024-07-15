##This script was written by Luke Archer (2024) and can be used to extract the gene ID column from a bedtools intersect output file using pandas 

##import pandas and os
import pandas as pd

import os

##Define user input
bedfile_path = input("Please enter the name of the bedtools intersect output file: ")

##print to screen to separate the file path from the filename
print(os.path.basename(bedfile_path).split("/")[-1])
##
bedfile_name=os.path.basename(bedfile_path)
infile=os.path.splitext(bedfile_name)

##The basename minus the suffix/ending is the first index of the input file
basename=infile[0]
suffix=infile[1]
#####

##Open an empty output file
output = []
##use pd.read_csv() to read the input file and use only the 12th column (index[11]) containing the gene IDs
infile = pd.read_csv(bedfile_path,sep="\t",usecols=[11])

#convert to pandas dataframe using pd.DataFrame
infile = pd.DataFrame(infile)
##use .to_csv() method to write a new output file as a .tsv file (sep="\t)
output=infile.to_csv(f"{basename}.IDs{suffix}",sep="\t")



