setwd("/Users/alexanderkissonergis/Documents/GitHub/ecologicalGenomics/transcriptomics/myresults/")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
BiocManager::install("vsn") 
BiocManager::install("hexbin")
### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

#############################
### Import our data
############################

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

#############################
### Explore Data Distributions
############################

# Let's see how many reads we have from each sample
colSums(countsTableRound)
mean(colSums(countsTableRound))

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,21000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] 2245.401
median(rowSums(countsTableRound)) # 117

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,1000), ylim=c(0,50000),breaks=10000)

#############################
### Start working with DESeq
############################

#### Create a DESeq object and define the experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

#############################
# CHECK THE QUALITY OF THE DATA BY SAMPLE CLUSTERING AND VARIATION
#############################
#install.packages("pheatmap")
#install.packages("RColorBrewer")

library("pheatmap")
library("vsn")

ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))


sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#AM_F11_Rep3 and OA_F2_Rep2 could be outliers (from resulting heat map and association above)

#####################################
# RUN PCA
####################################

# PCA to visualize global gene expression patterns

# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData,"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


####################################
# RUN A MORE ADVANCED PCA
####################################

# Let's plot the PCA by generation in four panels

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  

dataF0 <- subset(data, generation == 'F0')


F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  #ylim(-40, 25) + xlim(-50, 50)+ # new assembly limits
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  ##theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F0

#png("PCA_F0.png", res=300, height=5, width=5, units="in")

#ggarrange(F0, nrow = 1, ncol=1)

#dev.off()

################# F2

dataF2 <- subset(data, generation == 'F2')

F2 <- ggplot(dataF2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23), labels = c("Ambient", "Acidification","Warming"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A"), labels = c("Ambient", "Acidification","Warming"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F2

# png("PCA_F2.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F2, nrow = 1, ncol=1)
# 
# dev.off()

# Yes - F2 is missing one ambient replicate

################################ F4

dataF4 <- subset(data, generation == 'F4')

F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  #ylim(-20, 10) + xlim(-40, 25)+  # zoom with filtered assembly
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  #theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F4


# png("PCA_F4.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F4, nrow = 1, ncol=1)
# 
# dev.off()


################# F11

dataF11 <- subset(data, generation == 'F11')

F11 <- ggplot(dataF11, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-45, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,24), labels = c("Ambient", "OWA"))+
  scale_fill_manual(values=c('#6699CC', "#CC3333"), labels = c("Ambient", "OWA"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F11


# png("PCA_F11.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F11, nrow = 1, ncol=1)
# 
# dev.off()

ggarrange(F0, F2, F4, F11, nrow = 2, ncol=2)


#####################################
#Explore differential expression with full model, GE ~ generation + treatment
#####################################

resultsNames(dds)
#[1] "Intercept"            "generation_F11_vs_F0" "generation_F2_vs_F0"  "generation_F4_vs_F0" 
#[5] "treatment_OA_vs_AM"   "treatment_OW_vs_AM"   "treatment_OWA_vs_AM" 

resAM_OWA <- results(dds, name="treatment_OWA_vs_AM", alpha=0.05)

resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  

summary(resAM_OWA)
# log2 fold change (MLE): treatment OWA vs AM 
# Wald test p-value: treatment OWA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN22748_c0_g1::TRINITY_DN22748_c0_g1_i4::g.47585::m.47585   209.899      -1.799165  0.282811  -6.36172 1.99509e-10 2.75671e-06
# TRINITY_DN3600_c0_g1::TRINITY_DN3600_c0_g1_i2::g.16079::m.16079      91.998      -2.256764  0.355515  -6.34788 2.18302e-10 2.75671e-06
# TRINITY_DN505_c0_g1::TRINITY_DN505_c0_g1_i6::g.3715::m.3715         570.049      -0.698475  0.112387  -6.21493 5.13460e-10 4.32265e-06
# TRINITY_DN14069_c0_g1::TRINITY_DN14069_c0_g1_i2::g.39418::m.39418   193.197       3.555491  0.625276   5.68627 1.29844e-08 6.55866e-05
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733            699.038       2.606874  0.457072   5.70342 1.17424e-08 6.55866e-05
# TRINITY_DN8532_c1_g1::TRINITY_DN8532_c1_g1_i3::g.30152::m.30152     184.374       0.663567  0.119300   5.56217 2.66438e-08 1.12153e-04

# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 99, 0.39%
# LFC < 0 (down)     : 75, 0.3%
# outliers [1]       : 4, 0.016%
# low counts [2]     : 0, 0%
# (mean count < 17)

resAM_OW <- results(dds, name="treatment_OW_vs_AM", alpha=0.05)

head(resAM_OW)  

summary(resAM_OW)

# log2 fold change (MLE): treatment OW vs AM 
# Wald test p-value: treatment OW vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE         stat    pvalue      padj
# <numeric>      <numeric> <numeric>    <numeric> <numeric> <numeric>
#   TRINITY_DN1938_c0_g1::TRINITY_DN1938_c0_g1_i18::g.9836::m.9836     259.9095    1.08333e-01  0.129291  8.37898e-01  0.402088  0.999278
# TRINITY_DN5911_c0_g1::TRINITY_DN5911_c0_g1_i17::g.23209::m.23209   138.4247    1.33678e-02  0.108440  1.23275e-01  0.901890  0.999278
# TRINITY_DN209_c2_g1::TRINITY_DN209_c2_g1_i3::g.1533::m.1533         45.6569   -8.10373e-03  0.187617 -4.31930e-02  0.965548  0.999984
# TRINITY_DN10154_c0_g2::TRINITY_DN10154_c0_g2_i1::g.33548::m.33548   29.6715   -4.58357e-01  0.370513 -1.23709e+00  0.216055        NA
# TRINITY_DN6614_c0_g1::TRINITY_DN6614_c0_g1_i5::g.25280::m.25280     53.9376    2.36737e-05  0.600981  3.93917e-05  0.999969  0.999991
# TRINITY_DN6542_c3_g1::TRINITY_DN6542_c3_g1_i3::g.24952::m.24952     40.5470    5.04384e-02  0.249684  2.02009e-01  0.839910  0.999278
# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.036%
# LFC < 0 (down)     : 8, 0.032%
# outliers [1]       : 4, 0.016%
# low counts [2]     : 2449, 9.7%
# (mean count < 31)

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="TRINITY_DN1938_c0_g1::TRINITY_DN1938_c0_g1_i18::g.9836::m.9836", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
resAM_OW <- resAM_OW[order(resAM_OW$padj),]
head(resAM_OW)  

summary(resAM_OW)


##################################
#Run another model that focuses on treatment effects for each generation
##################################

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 67916     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_sub <- subset(dds, select = generation == 'F0')
dim(dds_sub)
# [1] 25260    12

# Perform DESeq2 analysis on the subset
dds_sub <- DESeq(dds_sub)

resultsNames(dds_sub)

res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

# log2 fold change (MLE): treatment OW vs AM 
# Wald test p-value: treatment OW vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130           32956.030        3.74319  0.292239   12.8086 1.46698e-37 3.63048e-33
# TRINITY_DN33_c0_g1::TRINITY_DN33_c0_g1_i18::g.235::m.235           4622.325        3.72565  0.296419   12.5689 3.13154e-36 3.87496e-32
# TRINITY_DN9365_c0_g1::TRINITY_DN9365_c0_g1_i1::g.32056::m.32056    3814.681        3.62002  0.292971   12.3562 4.50893e-35 3.71957e-31
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i18::g.1704::m.1704       4670.245        3.72840  0.317074   11.7588 6.36643e-32 3.20614e-28
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i13::g.1700::m.1700       5937.496        3.30257  0.280896   11.7573 6.47758e-32 3.20614e-28
# TRINITY_DN10209_c0_g2::TRINITY_DN10209_c0_g2_i3::g.33784::m.33784   279.313        3.00313  0.264004   11.3753 5.55141e-30 2.28977e-26

summary(res_F0_OWvAM)

# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3468, 14%
# LFC < 0 (down)     : 2049, 8.1%
# outliers [1]       : 22, 0.087%
# low counts [2]     : 490, 1.9%
# (mean count < 20)

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_sub, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7)
p

# We can make an MA plot
plotMA(res_F0_OWvAM, ylim=c(-6,6))

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By environment
vsd <- vst(dds_sub, blind=FALSE)

topgenes <- head(rownames(res_F0_OWvAM),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_sub)[,c("generation","treatment")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)

#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM

#################################################################

# For OW vs AM
res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)
res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)
res_F0_OWvAM <- res_F0_OWvAM[!is.na(res_F0_OWvAM$padj),]
degs_F0_OWvAM <- row.names(res_F0_OWvAM[res_F0_OWvAM$padj < 0.05,])

# For OA vs AM
res_F0_OAvAM <- results(dds_sub, name="treatment_OA_vs_AM", alpha=0.05)
res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)
res_F0_OAvAM <- res_F0_OAvAM[!is.na(res_F0_OAvAM$padj),]
degs_F0_OAvAM <- row.names(res_F0_OAvAM[res_F0_OAvAM$padj < 0.05,])

# For OWA vs AM
res_F0_OWAvAM <- results(dds_sub, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
#could further subset this based on the "stat" value to split by upregulated and downregulated (negative)
head(res_F0_OWAvAM)

summary(res_F0_OWAvAM)
res_F0_OWAvAM <- res_F0_OWAvAM[!is.na(res_F0_OWAvAM$padj),]
degs_F0_OWAvAM <- row.names(res_F0_OWAvAM[res_F0_OWAvAM$padj < 0.05,])

library(eulerr)

# Total
length(degs_F0_OAvAM)  # 602
length(degs_F0_OWvAM)  # 5517 
length(degs_F0_OWAvAM)  # 3918

# Intersections
length(intersect(degs_F0_OAvAM,degs_F0_OWvAM))  # 444
length(intersect(degs_F0_OAvAM,degs_F0_OWAvAM))  # 380
length(intersect(degs_F0_OWAvAM,degs_F0_OWvAM))  # 2743

intWA <- intersect(degs_F0_OAvAM,degs_F0_OWvAM)
length(intersect(degs_F0_OWAvAM,intWA)) # 338

# Number unique

602-444-380+338 # 116 OA
5517-444-2743+338 # 2668 OW 
3918-380-2743+338 # 1133 OWA

444-338 # 106 OA & OW
380-338 # 42 OA & OWA
2743-338 # 2405 OWA & OW


# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, "OW&OWA" = 2405, "OA&OW&OWA" = 338))


plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check
2668+2405+338+106 # 5517, total OW
1133+2405+338+42 # 3918, total OWA
116+42+106+338    # 602, total OA

####################################
# Work for 10.25.23
####################################

# Load the package
#install.packages("WGCNA")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("impute")

library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot) 
library(gridExtra)
library(Rmisc) 

#install.packages("remotes")
#remotes::install_github("kevinblighe/CorLevelPlot")

####################
#1. Import the counts matrix and metadata and filter using DESeq2
####################
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)
#67916, 38

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order

# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 6
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules.

cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs

head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module