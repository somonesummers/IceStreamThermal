%Plot Figure 3_1a
subplot(1,3,1)
cd figure3_1/
clear all; clc; path(pathdef);
addpath('..')
load data
plotBasalStressVerification(0,0,m,y,dz,Z,rho,g,alpha,tau_base,...
                            mu,u,MTP,'../data/schoofSoln_3_1a')

%Plot Figure 3_2a
subplot(1,3,2)
cd ../figure3_2/
clear all; clc; path(pathdef);
addpath('..')
load data
plotBasalStressVerification(0,0,m,y,dz,Z,rho,g,alpha,tau_base,...
                            mu,u,MTP,'../data/schoofSoln_3_2a')

%Plot Figure 3_3a
subplot(1,3,3)
cd ../figure3_3/
clear all; clc; path(pathdef);
addpath('..')
load data
plotBasalStressVerification(0,0,m,y,dz,Z,rho,g,alpha,tau_base,...
                            mu,u,MTP,'../data/schoofSoln_3_3a')
cd ..
setFontSize(13)