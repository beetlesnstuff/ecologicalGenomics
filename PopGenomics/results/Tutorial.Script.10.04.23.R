setwd("/Users/alexanderkissonergis/Documents/GitHub/ecologicalGenomics/PopGenomics/results/")

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

bio <- getData("worldclim",var="bio",res=10)

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)
