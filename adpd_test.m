
hc = xlsread('E:\BioSpyder\human_read_count.csv');
sum(sum(hc))/(38516728/4)*100%810886
mc = xlsread('E:\BioSpyder\mouse_read_count.csv');
sum(sum(mc))/(38516728/4)*100

hu = xlsread('E:\BioSpyder\human_umi_count.csv');
sum(sum(hu))/sum(sum(hc))*100
mu = xlsread('E:\BioSpyder\mouse_umi_count.csv');
sum(sum(mu))/sum(sum(mc))*100

b = xlsread('E:\BioSpyder\barcode_quantification.csv');
sum(b)/(38516728/4)*100

KLD = est_rel_entro_HJW;
% xlswrite('medNor_AD_KLD.csv',KLD');
 %xlswrite('medNor_PD_KLD.csv',KLD');
 %xlswrite('medNor_ADPD_KLD.csv',KLD');
%csvwrite('medNor_PD_KLD.csv',KLD');
%csvwrite('5samples_AD_KLD.csv',KLD');
%csvwrite('5samples_ADb_KLD.csv',KLD');
%csvwrite('5samplesAD_allCT_KLD.csv',KLD');
csvwrite('5samplesADc_KLD.csv',KLD');

P = xlsread('E:\BioSpyder\ADPD\medNor_PD_KLD.csv');
AP = xlsread('E:\BioSpyder\ADPD\medNor_ADPD_KLD.csv');
A = xlsread('E:\BioSpyder\ADPD\medNor_AD_KLD.csv');

A1=A';
find(A1>1.3)
figure,plot(A1,'*')
ylim([0 10])
max(A1)

P1 = P';
find(P1>1.2)
figure,plot(P1,'*')
ylim([0 10])
max(P1)

AP1=AP';
find(AP1>0.78)
figure,plot(AP1,'*')
ylim([0 10])
max(AP1)


%find(kld>2.5)

%2037        2324        3297


%pca(A)