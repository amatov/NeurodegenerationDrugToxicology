setwd ('C:\Users\Alex Matov\Documents\SeqData')
setwd("~/Bios2115")
library(FactoMineR)

A = read.csv('A_12.csv')
R = read.csv('R_12.csv')
T = read.csv('T_12.csv')

res.pca = PCA(A[,2:25])
res.pca = PCA(t(A[,2:25]))
res.pca = PCA(t(R[,2:25]))
res.pca = PCA(t(T[,2:25]))





BIOS2298_gene_counts

gene_counts <- read.csv("C:\Users\Alex Matov\OneDrive - biospyder.com\Documents\SeqData\BIOS2298_gene_counts.csv")
gene_counts <- read.csv("BIOS2298_gene_counts.csv")

gene_counts <- read.csv("C:\Users\Alex Matov\OneDrive - biospyder.com\Documents\SeqData\BIOS2298_gene_counts.csv")

#AP/PD
Rcounts <- read.csv("BIOS2180_gene_counts.csv")

library(FactoMineR)

library(Factoshiny)
library(missMDA)
library(FactoInvestigate)

library(Rcmdr)

data(decathlon)
res.pca = PCA(decathlon, quanti.sup=11:12,quali.sup=13)
library(Factoshiny)
resshiny = PCAshiny(res.pca)

res.pca = PCA(ds[1:3386,2:97])
plot(res.pca)
dst = t(ds)
res.pca = PCA(dst[2:97,])
