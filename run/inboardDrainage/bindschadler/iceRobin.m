clear all; clc; %path(pathdef);
addpath('../../../lib/','../../../lib/vis')
fID = fopen('log.txt','w');

%% %% %% %% %%     Initialization     %% %% %% %% %%

%%%%%     Runtime variables     %%%%%
%Domain
m = 1001;
n = 27; %m,n: dimensions of triangulation
Y = 32e3;
Z = 865; %Y,Z: dimensions of domain[m]
Z_till = 0; %Z_till: till thickness [m]
%Time
endT = 0; %endT: final time [s]
dt = 1; %dt: time step[s]
%Advection
v_0 = 0; %v_0: advection rate of ice from ridge [m/s]
a = 0; %a: accumulation rate [m/yr]
%Ice Rheology
iceRheol = 3; %iceRheol: ice rheology choice []
%1: Newtonian (1e14)        %2: Temp-Independent Glen's
%3: Temp-Dependent Glen's   %4: Temp-Dependent Goldsby
%Margin Solve Type
marginSolveType = 1;
%1: self-consistent margin  %2: choose margin location
%Coupling
omegaT = 0.3; %omegaT: thermal relaxation parameter []
omegaM = 0.5; %omegaR: rheology relaxation parameter []
tol = 1e-3; %tol: rheological and thermal error tolerance []
%Initialize constants
[rho,g,alpha,L,G_base,k1,k2,c1,c2,T_m,T_atm,...
 dy,dz,y,z,zT,nT,t,north,south,east,west,...
 northT,southT,eastT,westT,w,v] = initConstants(m,n,Y,Z,Z_till,endT,dt,v_0,a);
alpha = 0.0015;
G_base = 58.3e-3;

%%%%%     Pseudo-Initial Conditions     %%%%%
u(:,1) = zeros(m*n,1); %u: velocity[m/s]
T(:,1) = 260.15*ones(m*nT,1); %T: temperature[K]
mu(:,1) = 1e14*ones(m*n,1); %mu: effective viscosity[Pa s]
k(:,1) = k1*exp(-k2*1e-3.*T); %k: thermal conductivity[W/m*K]
c(:,1) = c1 + c2.*T; %c: specific heat[J/kg*K]

%%%%%     Basal Strength Evolution     %%%%%
drivingRamp = 1;
enforcedMargin = 16e3; %enforcedMargin: margin location [m]
tau_base = -175e3*exp(-abs(y-16.2e3)/0.25e3)...
           -50e3*exp(-abs(y-12e3)/0.15e3)...
           -56e3*exp(-abs(y-14e3)/0.1e3)...
           -50e3*exp(-abs(y-8e3)/0.15e3)...
           -45e3*exp(-abs(y-10e3)/0.1e3)...
           -0.52e3;


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

subplot(3,1,1)
plotSurfaceVelocity(m,n,y,u);
xlim([0 Y/1e3])
load ../digitizedVelocities.mat
hold on
plot((-profileI(:,1)+enforcedMargin)/1e3,profileI(:,2),'rx');

subplot(3,1,2)
plotBasalStress(0,0,m,y,dz,Z,rho,g,alpha,tau_base,mu,u,0)

subplot(3,1,3)
plotTemperatureField(m,nT,y,zT,T,T_m)

%%%   Plot Velocity Only   %%%
plot(y(1:m)/Y,u(m*n-m+1:m*n)/max(u(m*n-m+1:m*n)),'k','LineWidth',4)
hold on
plot((-profileI(:,1)+16.2e3)/Y,profileI(:,2)/(3.1557e7*max(u(m*n-m+1:m*n))),'r.','MarkerSize',20);
axis([0 0.52 0 1.1])
set(gca, 'xdir','reverse')
set(gcf,'color','w')

%%%   Plot Strength Only   %%%
h = plot(y(1:m)/Y,-(tau_base(1:m)./(Z*rho*g*sin(alpha))),'--','LineWidth',1.3);
set(h(1),'color',[0.7,0.7,0.7])
hold on
plot(y(1:m)/Y,-(mu(1:m)./(Z*rho*g*sin(alpha))).*(u(1:m)-u(m+1:2*m))/dz,'k','LineWidth',3)
axis([0 0.52 0 15])
set(gca, 'xdir','reverse')
set(gcf,'color','w')
