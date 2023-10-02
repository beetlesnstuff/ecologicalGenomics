setwd("/Users/alexanderkissonergis/Documents/GitHub/ecologicalGenomics/PopGenomics/results")

library(RcppCNPy)

list.files()

# read in the selection stats

s <- npyLoad("allRS_poly.selection.npy")

pval <- as.data.frame(1-pchisq(s,1))
head(pval)

names(pval) =c("p_PC1", "p_PC2")

#associate these with our SNP meta-data

p <- read.table("allRS_poly_mafs.sites", sep="\t", header=T, stringsAsFactors = T)

p_filtered = p[which(p$kept_sites==1),]

## make manhattan plot
plot(-log10(pval$p_PC1),
     col=p_filtered$chromo,
     xlab="Position",
     ylab="-log10(p-value)",
     main="Selection outliers: pcANGSD e=1 (K2)")

sel_contig <- p_filtered[which(pval==min(pval$p_PC1)), c("chromo", "position")]

cutoff=1e-3   # equals a 1 in 5,000 probability
outlier_contigs <- p_filtered[which(pval<cutoff),c("chromo","position")]
outlier_contigs

# how many outlier loci < the cutoff?
dim(outlier_contigs)[1]

# how many unique contigs harbor outlier loci?
length(unique(outlier_contigs$chromo))