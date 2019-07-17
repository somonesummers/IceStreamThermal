clear all; clc; path(pathdef);
addpath('../../../../lib/','../../../../lib/vis')
fID = fopen('log.txt','w');

%% %% %% %% %%     Initialization     %% %% %% %% %%

%%%%     Runtime variables     %%%%%
%Domain
m = 1001;
n = 11; %m,n: dimensions of triangulation
Y = 120e3;
Z = 700; %Y,Z: dimensions of domain[m]
Z_till = 0; %Z_till: till thickness [m]
%Time
endT = 0; %endT: final time [s]
dt = 5e7; %dt: time step[s]
%Advection
v_0 = 0; %v_0: advection rate of ice from ridge [m/s]
a = 1; %a: accumulation rate [m/yr]
%Ice Rheology
iceRheol = 3; %iceRheol: ice rheology choice []
%1: Newtonian (1e14)        %2: Temp-Independent Glen's
%3: Temp-Dependent Glen's   %4: Temp-Dependent Goldsby
%Margin Solve Type
marginSolveType = 1;
%1: self-consistent margin  %2: choose margin location
%Coupling
omegaT = 0.2; %omegaT: thermal relaxation parameter []
omegaM = 0.5; %omegaR: rheology relaxation parameter []
tol = 1e-3; %tol: rheological and thermal error tolerance []
%Initialize constants
[rho,g,alpha,L,G_base,k1,k2,c1,c2,T_m,T_atm,...
 dy,dz,y,z,zT,nT,t,north,south,east,west,...
 northT,southT,eastT,westT,w,v] = initConstants(m,n,Y,Z,Z_till,endT,dt,v_0,a);
alpha = 0.0007;
%%%%%%%%%%     Note: need 10*b in thermalModel.m     %%%%%%%%%%

%%%%%     Pseudo-Initial Conditions     %%%%%
u(:,1) = zeros(m*n,1); %u: velocity[m/s]
T(:,1) = 260.15*ones(m*nT,1); %T: temperature[K]
mu(:,1) = 1e14*ones(m*n,1); %mu: effective viscosity[Pa s]
k(:,1) = k1*exp(-k2*1e-3.*T); %k: thermal conductivity[W/m*K]
c(:,1) = c1 + c2.*T; %c: specific heat[J/kg*K]

%%%%%     Basal Strength Evolution     %%%%%
drivingRamp = 1*ones(1,length(t));
enforcedMargin = 30e3*ones(1,length(t)); %enforcedMargin: margin location [m]
for i = 1:length(t)
tau_base(:,i) = -40e3*exp((y-61.0e3)/1.2e3)...
                -(2.3e3+79e3*(t(i)/1e9))*exp(-abs(y-33.2e3)/(10.0e3-(t(i)/1e9)*8.3e3))...
                -2.0e3;
% tau_base(:,i) = -150e3*exp(-abs(y-61e3)/0.4e3)...
%                 -75e3*(t(i)/1e9)*exp(-abs(y-31.44e3)/0.4e3)...
%                 -3.07e3 + (t(i)/1e9)*1.6e3;
end


%% %% %% %% %%     Solution     %% %% %% %% %%

%%%%%     Steady-State Initial Conditions     %%%%%
fprintf(fID,'*****     Time: %d     *****\n',t(1));
[u(:,1),T(:,1),mu(:,1),k(:,1),c(:,1),MTP(:,1)] = ...
    marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,u(:,1),T(:,1),mu(:,1),...
                k(:,1),c(:,1),tau_base(:,1),enforcedMargin(:,1),north,...
                south,east,west,northT,southT,eastT,westT,iceRheol,...
                marginSolveType,rho,g,alpha,drivingRamp(1),k1,k2,c1,c2,...
                T_m,T_atm,G_base,w,v,omegaM,omegaT,tol,0);

%%%%%     Time Loop     %%%%%
for nt = 2:length(t)
    fprintf(fID,'*****     Time: %d     *****\n',t(nt));
    [u(:,nt),T(:,nt),mu(:,nt),k(:,nt),c(:,nt),MTP(:,nt)] = ...
        marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,u(:,nt-1),T(:,nt-1),...
                    mu(:,nt-1),k(:,nt-1),c(:,nt-1),tau_base(:,nt),...
                    enforcedMargin(:,nt),north,south,east,west,northT,...
                    southT,eastT,westT,iceRheol,marginSolveType,rho,g,...
                    alpha,drivingRamp(nt),k1,k2,c1,c2,T_m,T_atm,G_base,...
                    w,v,omegaM,omegaT,tol,1);
    plotTemperatureField(m,nT,y,zT,T(:,end),T_m);
    getframe;
end


%% %% %% %% %%     Visualization     %% %% %% %% %%
save data.mat
for i = 1:length(t)
    subplot(3,1,1)
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
    str = sprintf('%f years before present',t(i)/31622400);
    title(str)

    subplot(3,1,2)
    plotBasalStress(0,0,m,y,dz,Z,rho,g,alpha,tau_base(:,i),mu(:,i),u(:,i),MTP(:,i))
    axis([11.3 70 0 30])
    set(gca, 'xdir','reverse')

    subplot(3,1,3)
    plotTemperatureField(m,nT,y,zT,T(:,i),T_m)
    xlim([11.3 70])
    set(gca, 'xdir','reverse')
    setFontSize(16)
    set(gcf,'color','w');
    anim(i) = getframe(gcf);
end

%%%   Plot Surface Velocity Only   %%%
frame_1 = 1;
frame_2 = 5;
for i = frame_1:frame_2
    hold on
    %Find the crevassed region
    %if i == 21
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
    %end
    plot(y(1:m)/1000,3.1557e7*u(m*n-m+1:m*n,i),'k','LineWidth',4)
    set(gca,'Layer','top')
    axis([11.3 70 0 600])
    set(gca, 'XTick', [11.3,20:10:70])
    set(gca,'LineWidth',1.2)
    set(gca, 'xdir','reverse')
    set(gcf,'color','w')
    anim(i) = getframe;
    clf
end


%%%   Plot Strength Profile Only   %%%
for i = [6,13,17,21]
    hold on
    m2 = find(y==MTP(1));
    plot(y(1:m2)/1e3,-(tau_base(1:m2,i)./(Z*rho*g*sin(alpha))),...
         'Color',0.75*[1-(i-frame_1+1)/(frame_2-frame_1)...
         1-(i-frame_1+1)/(frame_2-frame_1)...
         1-(i-frame_1+1)/(frame_2-frame_1)],'LineWidth',2.2)
    axis([11.3 70 0 20])
    set(gca, 'XTick', [11.3,20:10:70])
    set(gca,'LineWidth',1.2)
    set(gca, 'xdir','reverse')
    set(gcf,'color','w')
    anim(i) = getframe;
end

%%%   Plot Strength & Stress Profiles   %%%
for i = 1:length(t)
    hold off
    h = plot(y(1:m)/1e3,-(tau_base(1:m,i)./(Z*rho*g*sin(alpha))),'--','LineWidth',1.3);
    set(h(1),'color',[0.7,0.7,0.7])
    hold on
    plot(y(1:m)/1e3,-(mu(1:m,i)./(Z*rho*g*sin(alpha))).*(u(1:m,i)-u(m+1:2*m,i))/dz,'k','LineWidth',3)
    axis([11.3 70 0 30])
    set(gca, 'XTick', [11.3,20:10:70])
    set(gca, 'xdir','reverse')
    set(gcf,'color','w')
    anim(i) = getframe;
end

%%%   Animate Strength Profile Only   %%%
frame_1 = 1;
frame_2 = length(t);
for i = frame_1:frame_2
    %hold on
    plot(y/1e3,-tau_base(:,i)/(rho*g*Z*alpha),'k','LineWidth',4)
    axis([11.3 70 0 32])
    set(gca, 'XTick', [11.3,20:10:70])
    set(gca, 'YTick', [0:10:30,32])
    set(gca,'LineWidth',1.2)
    set(gca, 'xdir','reverse')
    set(gcf,'color','w')
    anim(i) = getframe;
end
