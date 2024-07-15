
#this script can be used to output mean depth statistics for an F2.best.practice filtered VCF
import statistics 

##use statistics
input=open("110624.depth.per.site.ldepth","r").readlines()[1:]
#output=open("MeanDepthStatistics.txt","w")
#define the depth column
depths=[int(line.split('\t')[2]) for line in input]
#define mean depth
mean_depth=statistics.mean(depths)
#make an upper limit cut off for the depth
cut_off=int(round(1.6*mean_depth))

#output.write(mean_depth)

#output.write(cut_off)

print(mean_depth)

print(cut_off)

#input.close()
#output.close()

#####
##lastline

