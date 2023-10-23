setwd("~/ADPD")
setwd("Z:/kld")


# read ad pd and adpd files
A = read.csv('AD_cpm_KLD.csv')
A = read.csv('AD_2_KLD.csv')

P = read.csv('PD_cpm_KLD.csv')
AP = read.csv('PDAD_cpm_KLD.csv')

A1=A[,2]
ind1=which(A1>0)
A11=A1[ind1]
# plot histograms
hA <- hist( A11, breaks = 200, plot = FALSE) # 
plot(hA, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,2000)) # 

P1=P[,2]
ind2=which(P1>0)
P11=P1[ind2]
# plot histograms
hP <- hist( P11, breaks = 200, plot = FALSE) # 
plot(hP, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,2000)) # 

AP1=AP[,2]
ind3=which(AP1>0)
AP11=AP1[ind3]
# plot histograms
hAP <- hist( AP11, breaks = 200, plot = FALSE) # 
plot(hAP, col = rgb(1,0,0,1/10),xlim = c(0,2), ylim = c(0,2000)) # 

# look for the cutoff /unimodal?
library(FactoMineR)

#load from the ADPD set only the genes from the KLD200 list

#res.pca = PCA(ds[1:3386,2:97])
#plot(res.pca)
#dst = t(ds)
#res.pca = PCA(dst[2:97,])
