clear all
clc

load data
for i = 1:length(t)
%%%%%     Plot surface velocity and crevassed region     %%%%%
    subplot(2,1,1)
    hold off
    %Find the crevassed region
    indL = find(diff(gradient(3.1557e7*u(m*n-m+1:m*n,i))/dy <= -0.01)==1);
    indR = find(diff(gradient(3.1557e7*u(m*n-m+1:m*n,i))/dy <= -0.01)==-1);
    for j = 1:length(indL)
        h(j) = area([y(indL(j)) y(indR(j))]./1e3,[600 600]);
        if y(indR(j)) < 40e3
        set(h(j),'FaceColor',[1,0.67,0.65],'EdgeColor',[1,0.67,0.65])
        else
        set(h(j),'FaceColor',[0.72,0.80,1],'EdgeColor',[0.72,0.80,1])
        end
        hold on
    end
    plotSurfaceVelocity(m,n,y,u(:,i));
    set(gca,'Layer','top')
    axis([11.3 70 0 600])
    set(gca,'YTick',[0,210,500])
    set(gca, 'xdir','reverse')

%%%%%     Plot basal stress and strength     %%%%%
    subplot(2,1,2)
    plotBasalStress(0,0,m,y,dz,Z,rho,g,alpha,tau_base(:,i),mu(:,i),u(:,i),MTP(:,i))
    axis([11.3 70 0 30])
    set(gca, 'xdir','reverse')
    setFontSize(16)
    set(gcf,'color','w');
    anim(i) = getframe(gcf);
    
%%%%%     Save animation as AVI     %%%%%
movie2avi(anim,'kambMigrationAnimation.avi','FPS',3)
end