## Load the package
library("FactoMineR")
library(xlsx)

setwd("E:/BioSpyder/Bios2115") # change accordingly

data <- read.xlsx("BIOS2115_bios2115_MCF7_ART_gene_countsLABELLED.xlsm")

data1 <- read.csv("BIOS2115_bios2115_plate4_gene_countsLABELLED.csv", sheet = 2)

data2 <- read.csv("test.txt", sheet=3)

d2 = data.frame(data2)


res.pca <- PCA(matrix)   
plot(res.pca)				       
