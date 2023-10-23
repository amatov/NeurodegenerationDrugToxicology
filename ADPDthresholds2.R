setwd("~/ADPD")
library("RColorBrewer")
#library(pvclust)
library(tsne)
tsne <- Rtsne(heatPD, dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
library(FactoMineR)

DAD = read.csv('Alzheimer.csv_cpm.csv', header = F)
DAD5 = DAD[,1:20]
DPD = read.csv('Parkinson.csv_cpm.csv', header = F)
DCT = read.csv('Control.csv_cpm.csv', header = F)
DCT5 = DCT[,1:20]
A = read.csv('AD_cpm_KLD.csv', header = F)
A1=A[,2]
indKLD_ADshort = which(A1>0.9113598)
P = read.csv('PD_cpm_KLD.csv', header = F)
P1=P[,2]
indKLD_PD = which(P1>1.3)#0.9113598)

KLD5_AD = read.csv('5samples_AD_KLD.csv', header = F)
indA5 = which(KLD5_AD>1)
# 2490  3028  5624  7869 10527 16605 16964 18074 19329 21136 21442
KLD5_ADb = read.csv('5samples_ADb_KLD.csv', header = F)
indA5 = which(KLD5_ADb>1)
KLD5_ADall = read.csv('5samplesAD_allCT_KLD.csv', header = F)
indA5 = which(KLD5_ADall>0.25)

KLD5_ADc = read.csv('5samplesADc_KLD.csv', header = F)
indA5 = which(KLD5_ADc>0.1)


DAD5 = DAD[indA5,1:20]
colnames(DAD5)<-c(1:20)
DCT5 = DCT[indA5,29:48]
DCTa = DCT[indA5,1:88]
#colnames(auxP)<-c('CT')
test10 = cbind(DAD5,DCTa)

test10 = cbind(DCT5, DAD5)
df<-scale(t(test10))
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col)



res.pca = PCA(t(test10))
# add test t-SNE
# add test UMAP


res.pca = PCA(DAD)
plot(res.pca)
dst = t(ds)
res.pca = PCA(dst[2:97,])
res.pca = PCA(DAD[indKLD_PD,])
res.pca = PCA(auxC)
res.pca = PCA(auxP)
res.pca = PCA(t(heatPD))

indKLD_PA = which(PA1>0.7212865)
auxP <- PAD[indKLD_PD,]
colnames(auxP)<-c(101:222)

##################
auxC <- DCT[indKLD_ADshort,]
dim(auxC)
#colnames(auxC)<-c(1:88)
hv2 <-matrix(,nrow=11,ncol=88)
hv2<-auxC
names(hv2) <- NULL

auxA = DAD[indKLD_ADshort,]
dim(auxA)

heatAD = cbind(auxA,hv2)

df<-scale(t(heatAD))
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col)
#########################

heatPD = cbind(auxP,auxC)

heatPD = cbind(PAD[indKLD_PD,],CAD[indKLD_PD,])
#hierarchical clusteirng on KLD200
###########################################################
# Hierarchical cluster analysis on a set of dissimilarities and methods for analyzing it. 
aux <- t(mafsQ2[,1:5])
hc <- hclust(aux)               
plot(hc) 

# A heatmap is another way to visualize hierarchical clustering; in the colored image, data values are transformed to color scale.
df<-scale(t(heatAD))
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
heatmap(df, scale = "none", col =  col)

# Hierarchical Clustering with P-Values via Multiscale BootstrapResampling
# https://cran.r-project.org/web/packages/pvclust/pvclust.pdf
result <- pvclust(t(mafsP2[,1:5000]), method.dist="cor", method.hclust="average", nboot=1000, parallel=TRUE)
plot(result)


#tSNE clustering on full gene set and on KLD200



# read ad pd and adpd and pdad CPM files
A = read.csv('AD_cpm_KLD.csv', header = F)
A1=A[,2]
write.csv(x=A1,'AD_cpm_KLD_v.csv')

A11= sort(A1)
ind1 = length(A11)
A111= A11[(ind1-199):ind1]
indKLD_AD = which(A1>0.5821707)
length(A1[indKLD_AD])
write.csv(x=indKLD_AD,'AD_cpm_KLD200.csv')

P = read.csv('PD_cpm_KLD.csv', header = F)
P1=P[,2]
write.csv(x=P1,'PD_cpm_KLD_v.csv')

P11= sort(P1)
ind2 = length(P11)
P111= P11[(ind2-199):ind2]
indKLD_PD = which(P1>0.9113598)
length(P1[indKLD_PD])

write.csv(x=indKLD_PD,'PD_cpm_KLD200.csv')

AP = read.csv('ADPD_cpm_KLD.csv', header = F)
AP1=AP[,2]
write.csv(x=AP1,'ADPD_cpm_KLD_v.csv')

AP11= sort(AP1)
ind3 = length(AP11)
AP111= AP11[(ind3-199):ind3]
indKLD_AP = which(AP1>0.6548049)
length(A1[indKLD_AP])

write.csv(x=indKLD_AP,'ADPD_cpm_KLD200.csv')

PA = read.csv('PDAD_cpm_KLD.csv', header = F)
PA1=PA[,2]
write.csv(x=PA1,'PDAD_cpm_KLD_v.csv')

PA11= sort(PA1)
ind4 = length(PA11)
PA111= PA11[(ind4-199):ind4]

indKLD_PA = which(PA1>0.7212865)
length(PA1[indKLD_PA])
write.csv(x=indKLD_PA,'PDAD_cpm_KLD200.csv')

# read ad pd and adpd files
A = read.csv('AD_KLD.csv')
P = read.csv('PD_KLD.csv')
AP = read.csv('ADPD_KLD.csv')

APN = read.csv('medNor_ADPD_KLD.csv')
PN = read.csv('medNor_PD_KLD.csv')
AN = read.csv('medNor_AD_KLD.csv')

#A1=A[,2]
#ind1=which(A1>0)
ind1=which(t(AN)>0)
A11=t(AN)[ind1]
#A11=A1[ind1]
# plot histograms
hA <- hist( A11, breaks = 200, plot = FALSE) # 
plot(hA, col = rgb(1,0,0,1/10),xlim = c(0,10), ylim = c(0,1500)) # 

#P1=P[,2]
#ind2=which(P1>0)
#P11=P1[ind2]
ind2=which(t(PN)>0)
P11=t(PN)[ind2]
# plot histograms
hP <- hist( P11, breaks = 200, plot = FALSE) # 
plot(hP, col = rgb(1,0,0,1/10),xlim = c(0,10), ylim = c(0,1500)) # 

#AP1=AP[,2]
ind3=which(t(APN)>0)
AP11=t(APN)[ind3]
# plot histograms
hAP <- hist( AP11, breaks = 200, plot = FALSE) # 
plot(hAP, col = rgb(1,0,0,1/10),xlim = c(0,10), ylim = c(0,1500)) # 
#plot(hAP, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,2000)) # 

# look for the cutoff /unimodal?

# hierarchical clustering
hcTop20 <- rbind(cv2 , hv2)
df<-scale(hcTop20)
col <- colorRampPalette(brewer.pal(11, "RdYlBu"))(256)
hm <- heatmap(df, scale = "none", col =  col) 