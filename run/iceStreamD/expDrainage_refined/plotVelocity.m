clear all
clc

set(0,'defaultLineLineWidth',4)
set(0,'DefaultLineMarkerSize',10)
% subplot(2,1,1)
for i = 6:-1:2
    col = gray(8);
    col = [col(:,1),zeros(8,1),zeros(8,1)];
    col = col(2:8,:);
    cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
    offset = [15;15.4;17.8;19.6;24.3;30.1];
    load ../profileVelData.mat
    hold on
    cd(cases(i,:))
    load data
    evaluation = sprintf('plot(-%s(:,1)+%d,%s(:,2),''.'',''color'',col(%d,:))',cases(i,:),offset(i,:)-MTP/1e3,cases(i,:),i);
    eval(evaluation)
    cd ..
end
box on
axis([0 37 0 700])
setFontSize(12)
ylabel('Surface Velocity [m/yr]')
xlabel('Distance from Stream Center [km]')

for i = 6:-1:2
    col = gray(8);
    col = [col(:,1),zeros(8,1),zeros(8,1)];
    col = col(2:8,:);
    cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
    cd(cases(i,:))
    load data
    plot(y(1:m)/1000-MTP/1e3,3.1557e7*u(m*n-m+1:m*n),'Color',col(i,:))
    hold on
    cd ..
    clear all
end
set(gcf,'color','w');
box on
axis([-37 0 0 700])
setFontSize(12)
ylabel('Surface Velocity [m/yr]')
xlabel('Distance from Stream Center [km]')
set(gca,'linewidth',2,'Layer','top')

% subplot(2,1,2)
% for i = 6:-1:2
%     col = gray(8);
%     col = col(1:6,:);
%     cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
%     offset = [15;15;17.5;20;25;30];
%     load ../profileVelData.mat
%     hold on
%     evaluation = sprintf('plot(-%s(:,1)+%d,%s(:,2),''.'',''color'',col(%d,:))',cases(i,:),offset(i,:),cases(i,:),i);
%     eval(evaluation)
% end
% box on
% axis([0 30 0 700])
% setFontSize(12)
% ylabel('Surface Velocity [m/yr]')
% xlabel('Distance from Stream Center [km]')
% 
% for i = 6:-1:2
%     col = gray(8);
%     cd ../frozenTemperate/
%     cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
%     cd(cases(i,:))
%     load data
%     plot(y(1:m)/1000,3.1557e7*u(m*n-m+1:m*n),'Color',col(i,:))
%     hold on
%     cd ..
%     clear all
% end
% set(gcf,'color','w');
% box on
% axis([0 30 0 700])
% setFontSize(12)
% ylabel('Surface Velocity [m/yr]')
% xlabel('Distance from Stream Center [km]')
% set(gca,'linewidth',2,'Layer','top')
% 
% cd ../drainage2
set(0,'defaultLineLineWidth',1)
set(0,'DefaultLineMarkerSize',1)