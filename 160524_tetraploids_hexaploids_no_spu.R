#Script for conducting PCA on mixed ploidy VCF files
#Tuomas Hämälä, April 2024, Modified by Luke Archer (May 2024).

setwd("/Users/lukearcher/Desktop/Individual_Project_files/VCF_files/")
library(vcfR)
library(ggplot2)
library(ggrepel)
library(gridExtra) #for making multipanel plots

##load the different VCFs - firstly the ld pruned VCF with all UK tetraploids &&  
##all UK hexaploids (INCLUDING SPU)
vcf <- read.vcfR("160524_ld_pruned_20PCTmis_maf005_allUKtets_allUKhex.vcf.gz")#Filtered and LD-pruned VCF file

##now load the second VCF with all UK tetraploids and all hexaploids except SPU
vcf_2 <-read.vcfR("160524_ld_pruned_20PCTmis_maf005_NO_SPU.vcf.gz")

#Transform the first VCF to numeric genotypes
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(y)as.numeric(as.character(y))))

#Transform the second VCF (without SPU) to numeric genotypes
df_2 <- extract.gt(vcf_2)
df_2[df_2 == "0|0"] <- 0
df_2[df_2 == "0|1"] <- 1
df_2[df_2 == "1|0"] <- 1
df_2[df_2 == "1|1"] <- 2
df_2[df_2 == "0/0"] <- 0
df_2[df_2 == "0/1"] <- 1
df_2[df_2 == "1/0"] <- 1
df_2[df_2 == "1/1"] <- 2
df_2[df_2 == "0/0/0/0"] <- 0
df_2[df_2 == "0/0/0/1"] <- 1
df_2[df_2 == "0/0/1/1"] <- 2
df_2[df_2 == "0/1/1/1"] <- 3
df_2[df_2 == "1/1/1/1"] <- 4
df_2[df_2 == "0/0/0/0/0/0"] <- 0
df_2[df_2 == "0/0/0/0/0/1"] <- 1
df_2[df_2 == "0/0/0/0/1/1"] <- 2
df_2[df_2 == "0/0/0/1/1/1"] <- 3
df_2[df_2 == "0/0/1/1/1/1"] <- 4
df_2[df_2 == "0/1/1/1/1/1"] <- 5
df_2[df_2 == "1/1/1/1/1/1"] <- 6
df_2<- data.frame(apply(df_2,2,function(y)as.numeric(as.character(y))))


############
##MODIFIED BY LUKE ARCHER ON 16/05/2024 
#REMOVE SAMPLES WITH > 20% MISSING DATA
mis <- apply(df,2,function(y)sum(is.na(y))/length(y))
df <- df[,mis <= 0.2]

mis_2 <- apply(df_2,2,function(y)sum(is.na(y))/length(y))
df_2 <- df_2[,mis <= 0.2]

#Calculate allele frequencies
x <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(y)sum(y,na.rm=T)/sum(x[!is.na(y)]))

#Calculate allele frequencies for 2nd vcf
x2 <- apply(df_2,2,max,na.rm=T)
p2 <- apply(df_2,1,function(y)sum(y,na.rm=T)/sum(x[!is.na(y)]))

#Removing individuals can change allele frequencies, so we make sure that maf >= 0.05
df <- df[p >= 0.05 & p <= 0.95,]
p <- p[p >= 0.05 & p <= 0.95]

#Removing individuals can change allele frequencies, so we make sure that maf >= 0.05
##do the same for the second VCF
df_2 <- df_2[p2 >= 0.05 & p2 <= 0.95,]
p2 <- p2[p2 >= 0.05 & p2 <= 0.95]

#Estimate a covariance matrix
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:i){
    cov[i,j] <- mean((df[,i]/x[i]-p)*(df[,j]/x[j]-p)/(p*(1-p)),na.rm=T)
    cov[j,i] <- cov[i,j]
  }	
}

#Estimate a covariance matrix for the NO SPU vcf
n2 <- ncol(df_2)
cov2 <- matrix(nrow=n2,ncol=n2)
for(i in 1:n2){
  for(j in 1:i){
    cov2[i,j] <- mean((df_2[,i]/x[i]-p2)*(df_2[,j]/x[j]-p2)/(p2*(1-p2)),na.rm=T)
    cov2[j,i] <- cov2[i,j]
  }	
}

#Do PCA on the matrix and plot
pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=x)

#Do PCA on the matrix and plot the NO SPU vcf
pc2 <- prcomp(cov2,scale=T)
xlab2 <- paste0("PC1 (",round(summary(pc2)$importance[2]*100),"%)")
ylab2 <- paste0("PC2 (",round(summary(pc2)$importance[5]*100),"%)")
pcs2 <- data.frame(PC1=pc2$x[,1],PC2=pc2$x[,2],id=colnames(df_2),ploidy=x2)

############
##MODIFIED SLIGHTLY TO ADD AN ELLIPSE TO CIRCLE THE SAMPLES BY PLOIDY
with_SPU<-ggplot(pcs, aes(PC1, PC2, color=as.factor(ploidy)))+
  geom_point(size=7)+
  stat_ellipse()+
  labs(x=xlab, y=ylab, color="Ploidy")+
  geom_text_repel(aes(label=id), size=2.5, force=30, color="black")+
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks=element_line(color="black",size=0.5),
        axis.title=element_text(size=12, color="black"),
        plot.title=element_text(size=14, color="black", hjust = 0.5),
        legend.text=element_text(size=11, color="black"),
        legend.title=element_text(size=12, color="black"),
        legend.key=element_blank(),
        aspect.ratio=1)
#############


################
##WITHOUT SPU
without_spu<-ggplot(pcs2, aes(PC1, PC2, color=as.factor(ploidy)))+
  geom_point(size=7)+
  stat_ellipse()+
  labs(x=xlab, y=ylab, color="Ploidy")+
  geom_text_repel(aes(label=id), size=3, force=30, color="black")+
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(color="black",linewidth=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks=element_line(color="black",size=0.5),
        axis.title=element_text(size=12, color="black"),
        plot.title=element_text(size=14, color="black", hjust = 0.5),
        legend.text=element_text(size=11, color="black"),
        legend.title=element_text(size=12, color="black"),
        legend.key=element_blank(),
        aspect.ratio=1)

#######
##make a multipanel plot with grid.arrange()
grid.arrange(with_SPU, without_spu,nrow=1)
