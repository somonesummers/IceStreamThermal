clear all
clc

set(0,'defaultLineLineWidth',4)
hold on
box on
for i = 2:6
    col = gray(8);
    col = col(1:6,:);
    cases = ['profileF';'profileE';'profileD';'profileC';'profileB';'profileA'];
    cd(cases(i,:))
    load data
    plot(MTP/1e3-y(y<MTP)/1e3,(i-1)*11-(tau_base(y<MTP)./(Z*rho*g*sin(alpha))),'Color',col(i,:))
    plot(0:35:35,(i-1)*11*ones(2,1),'k','LineWidth',1)
    cd ..
    clear all
    xlabel('Distance from Stream Center [km]')
    set(gca,'ytick',[])
    axis([0 30 11 66])
    set(gca, 'xdir','reverse')
    set(gcf,'color','w');
    setFontSize(12)
    set(gca,'linewidth',2,'Layer','top')
end
set(0,'defaultLineLineWidth',1)
