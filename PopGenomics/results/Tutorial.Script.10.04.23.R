setwd("/Users/alexanderkissonergis/Documents/GitHub/ecologicalGenomics/PopGenomics/results/")

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])

clim <- extract(bio,points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,clim)
str(angsd_coords_clim)

# Make the climate PCA:

clim_PCA = PCA(angsd_coords_clim[,15:33], graph=T)

# Get a screeplot of cliamte PCA eigenvalues

fviz_eig(clim_PCA)

# What is the climate PCA space our red spruce pops occupy?

fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

# Which variables show the strongest correlation on the first 2 climate PC axes?

dimdesc(clim_PCA)[1:2]

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio12"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)


# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)