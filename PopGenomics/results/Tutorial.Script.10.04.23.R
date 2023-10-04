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