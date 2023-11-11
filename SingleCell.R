r = 38516728 # number of all reads wc -l BIOS2319-S01_S1_L001_R1_001.fastq
hc = read.csv('human_read_count.csv')#, header = F) # 22537x116=2614292 hc[1:22537,2:116]
mc = read.csv('mouse_read_count.csv')#, header = F)
hu = read.csv('human_umi_count.csv')#, header = F)
mu = read.csv('mouse_umi_count.csv')#, header = F)
b = read.csv('barcode_quantification.csv')#, header = F)

#Rate of mapped human reads = sum (human_read_count.csv)/(READS/4)
dimH = dim(hc)      
shc = sum(unlist(hc[,2:dimH[2]])) # 810886
hmr = shc/(r/4)*100#8%

#Rate of mapped mouse reads = sum (mouse_read_count.csv)/(READS/4)
dimM = dim(mc)
smc = sum(unlist(mc[,2:dimM[2]]))#1340945
mmr = smc/(r/4)*100#14%

#Percent UMI human = sum (human_umi_count.csv)/sum of mapped human reads
dimHU = dim(hu)
shu = sum(unlist(hu[,2:dimHU[2]]))#233833
puh = shu/shc*100#29%

#Percent UMI mouse = sum (mouse_umi_count.csv)/sum of mapped mouse reads
dimMU = dim(mu)
smu =  sum(unlist(mu[,2:dimMU[2]]))#363913
pum = smu/smc*100#27%

#Barcoding rate = sum (barcode_quantification.csv)/(READS/4)
sb = sum(b[,2]) #8077328 # length barcoding vector 60881
br = sb/(r/4)*100#84%

#Number of found human cells = #columns (human_read_count.csv)
aux1 = dim(hc)
nhc = aux1[2]-1#115

#Number of found human cells = #columns (human_umi_count.csv)
aux11 = dim(hu)
nhc1 = aux11[2]-1#115

#Number of found mouse cells = #columns(mouse_read_count.csv)
aux2 = dim(mc)
nmc = aux2[2]-1#115

#Number of found mouse cells = #columns(mouse_umi_count.csv)
aux22 = dim(mu)
nmc1 = aux22[2]-1#115