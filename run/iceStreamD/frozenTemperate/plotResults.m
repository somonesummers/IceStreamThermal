clear all
clc

set(0,'defaultLineLineWidth',4)
set(0,'DefaultLineMarkerSize',10)
for i = 1:6
    col = lines(12);
    col = col(7:end,:);
    cases = ['profileA';'profileB';'profileC';'profileD';'profileE';'profileF'];
    offset = [29.5;25;20;17.5;15.5;15.5];
    load ../profileVelData.mat
    hold on
    evaluation = sprintf('plot(-%s(:,1)+%d,%s(:,2),''.'',''color'',col(%d,:))',cases(i,:),offset(i,:),cases(i,:),i);
    eval(evaluation)
end
for i = 1:6
    cases = ['profileA';'profileB';'profileC';'profileD';'profileE';'profileF'];
    cd(cases(i,:))
    load data
    plotSurfaceVelocity(m,n,y,u);
    cd ..
    clear all
end
set(gcf,'color','w');
box on
axis([0 30 0 700])
setFontSize(12)
set(gca,'linewidth',2,'Layer','top')