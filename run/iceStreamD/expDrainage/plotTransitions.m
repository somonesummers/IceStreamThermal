clear all
clc

set(0,'defaultLineLineWidth',4)
set(0,'DefaultLineMarkerSize',10)
center = [30,25,20,17.5,15,15];
for i = 1:6
    cases = ['profileA';'profileB';'profileC';'profileD';'profileE';'profileF'];
    cd(cases(i,:))
    load data
    MTPindex = find(u(1:m)*3.15e7 < 1e-10,1);
    TTPindex = find(T(1:m)-T_m < -1e-10,1);
    spacing(i) = y(TTPindex) - y(MTPindex);
    if spacing(i) < 0
        spacing(i) = 0;
    end
    subplot(6,1,i)
    hold on
    plotTemperatureField(m,nT,y,zT,T,T_m);
    axis([center(i) - 2.5, center(i) + 2.5, 0, Z])
    plot(y(TTPindex)/1e3,0,'kd','LineWidth',1)
    plot(y(MTPindex)/1e3,0,'ko','LineWidth',1)
    cd ..
end
figure
plot(1:6,spacing,'k--.','LineWidth',2,'MarkerSize',50);
xlabel('Profile')
ylabel('Spacing between transitions [m]')
axis([1 6 0 500])
set(gcf,'color','w');
box on
setFontSize(15)
set(gca,'linewidth',2,'Layer','top')