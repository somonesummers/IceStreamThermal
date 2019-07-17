function[] = plotBasalStressBenchmark(offset,eps,m,y,dz,Z,rho,g,alpha,...
                                      tau_base,mu,u,marginPosition,dataPath)
%Plot basal stress, smearing region, and basal strength profile normalized
%by driving stress

h = area([offset offset+eps]./1000,[30 30]);
set(h(1),'FaceColor',[0.9,0.95,1],'EdgeColor',[0.9,0.95,1])
hold on
h = plot(y(1:m)/1000,-(tau_base(1:m)./(Z*rho*g*sin(alpha))),'--',...
         y(1:m)/1000,ones(m,1),'--');
set(h(1),'color',[0.7,0.7,0.7])
set(h(2),'color',[0.7,0.7,0.7])
str = sprintf('load %s',dataPath);
eval(str);
plot(schoofSoln(:,1),schoofSoln(:,2),'r.','MarkerSize',20)
plot(y(1:m)/1000,-(mu(1:m)./(Z*rho*g*sin(alpha))).*(u(1:m)-u(m+1:2*m))/dz,...
     'k','LineWidth',2)
axis([0 max(y/1e3) 0 30])
xlabel('Distance from Stream Center [km]')
ylabel('Normalized Basal Stress')
str = sprintf('Margin Position = %1.2f km',marginPosition/1e3);
title(str)
set(gca,'Layer','top')
hold off

end