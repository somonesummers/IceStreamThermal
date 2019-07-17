%Plot Figure 4_1
subplot(3,2,1)
cd figure4_1/
clear all; clc; path(pathdef);
addpath('../../../lib/','../../../lib/vis')
load data
plotTemperatureField(m,nT,y,zT,T,T_m)
subplot(3,2,2)
hold on
load('../upBSurfVel.mat')
plot(y_ob,u_ob,'.r','MarkerSize',20)
plotSurfaceVelocity(m,n,y,u);
axis([0 20 0 500]), box on

%Plot Figure 4_2
subplot(3,2,3)
cd ../figure4_2/
clear all; clc; path(pathdef);
addpath('../../../lib/','../../../lib/vis')
load data
plotTemperatureField(m,nT,y,zT,T,T_m)
subplot(3,2,4)
hold on
load('../upBSurfVel.mat')
plot(y_ob+2.8,u_ob,'.r','MarkerSize',20)
plotSurfaceVelocity(m,n,y,u);
axis([0 20 0 500]), box on

%Plot Figure 4_3
subplot(3,2,5)
cd ../figure4_3/
clear all; clc; path(pathdef);
addpath('../../../lib/','../../../lib/vis')
load data
plotTemperatureField(m,nT,y,zT,T,T_m)
subplot(3,2,6)
hold on
load('../upBSurfVel.mat')
plot(y_ob+4.1,u_ob,'.r','MarkerSize',20)
plotSurfaceVelocity(m,n,y,u);
axis([0 20 0 500]), box on

cd ..
setFontSize(13)