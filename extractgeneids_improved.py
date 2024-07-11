##this script can be used to extract the gene ID column from a bedtools intersect output file
import pandas as pd

import os

##define user input
bedfile_path = input("Please enter the name of the bedtools intersect output file: ")

print(os.path.basename(bedfile_path).split("/")[-1])

bedfile_name=os.path.basename(bedfile_path)
infile=os.path.splitext(bedfile_name)

print(infile)

print(infile[0] + infile[1])

basename=infile[0]
suffix=infile[1]
#######
#use argparse() for the input and output arguments
##write some code that will go through a directory/folder of .tsv files and take the "basename" of the file and then correctly produce the output with the right basename
##in one of the files, the index[11] was wrong and the genee ID was not actually in this column, so find a way to extract the column regardless of position

##open an empty output file
output = []
##use pd.read_csv() to read the input file and use only the 12th column (index[11])
infile = pd.read_csv(bedfile_path,sep="\t",usecols=[11]) ##have to manually change this for ang_pyr_dan_5010.geneoverlaps.tsv because index[11] is not the right column

#convert to pandas dataframe
infile = pd.DataFrame(infile)
#print to screen to check the right column has been selected
##use .to_csv() method to produce a new output file
print(infile)
output=infile.to_csv(f"{basename}.IDs{suffix}",sep="\t")

#geneid=infile[11]

#output=output.write(geneid)

print(output)

