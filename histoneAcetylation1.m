function  nuc = histoneAcetylation1

% mirs = xlsread('C:\Users\dsa\Desktop\mirSet02Raw.xlsx');
nuc = imread('D:\LESHA\MPP8-Katoshka-Test2\araC\Katoshka__C_005_r_0005_c_0005_t_00000000_z_0000-00000000.tif');
%I=nuc(1:170,70:200);
I=nuc;
figure,imshow(I,[])

%points = orb(I);
points = detectSURFFeatures(I);
[features, valid_points] = extractFeatures(I, points);

figure
imshow(I,[])
hold on
plot(valid_points.selectStrongest(10),'showOrientation',true);
%plot(points,'ShowScale',false)
hold off



% figure,plot(mirs(1,:))
% hold on,plot(mirs(2,:),'red')
% hold on,plot(mirs(3,:),'green')
% hold on,plot(mirs(4,:),'yellow')


figure,plot(mirs(1,1:15))
hold on,plot(mirs(2,1:15),'red')
hold on,plot(mirs(3,1:15),'green')
hold on,plot(mirs(4,1:15),'yellow')