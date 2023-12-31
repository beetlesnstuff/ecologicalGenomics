setwd("/Users/alexanderkissonergis/Documents/GitHub/ecologicalGenomics/PopGenomics/results/")

list.files()

# First let's read in the diversity stats
theta <- read.table("2103_.thetas",sep="\t",header=T)

head(theta)

theta$tWsite = theta$tW/theta$nSites #scales the theta-W by the number of sites
theta$tPsite = theta$tP/theta$nSites #scales the theta-Pi by the number of sites

summary(theta)

# You can order the contig list to show you the contigs with the highest values of Tajima's D, or the lowest

head(theta[order(theta$Tajima, decreasing = TRUE),]) # top 10 Tajima's D values

head(theta[order(theta$Tajima, decreasing = FALSE),]) # bottom 10 Tajima's D values

#You can also look for contigs that have combinations of high Tajima's D and low diversity -- these could represent outliers for selection
#theta[which(theta$Tajima>1.5 & theta$tPsite<0.001),]

sfs<-scan('2103_.sfs')
sfs<-sfs[-c(1,which(sfs==0))]
sfs<-sfs/sum(sfs)

# Be sure to replace "9999" with your pop code in the "main" legend below
barplot(sfs,xlab="Chromosomes",
        names=1:length(sfs),
        ylab="Proportions",
        main="Pop 2103 Site Frequency Spectrum",
        col='blue')

# Put the nucleotide diversities, Tajima's D, and SFS into a 4-panel figure
par(mfrow=c(2,2))
hist(theta$tWsite, xlab="theta-W", main="Watterson's theta")
hist(theta$tPsite, xlab="theta-Pi", main="Pairwise Nucleotide Diversity")
hist(theta$Tajima, xlab="D", main="Tajima's D")
barplot(sfs,names=1:length(sfs),main='Site Frequency Spectrum')

sum(theta$nSites)



# To reset the panel plotting, execute the line below:
dev.off()