clear all; clc; path(pathdef);
addpath('../../lib/','../../lib/vis')
fID = fopen('log.txt','w');

%% %% %% %% %%     Initialization     %% %% %% %% %%
for kk = 1:20
%%%%     Runtime variables     %%%%%
%Domain
m = 1001;
n = 11; %m,n: dimensions of triangulation
Y = 100e3;
Z = 1e3; %Y,Z: dimensions of domain[m]
Z_till = 0; %Z_till: till thickness [m]
%Time
endT = 0; %endT: final time [s]
dt = 1; %dt: time step[s]
%Advection
v_0 = 0; %v_0: advection rate of ice from ridge [m/s]
a = 0; %a: accumulation rate [m/yr]
%Ice Rheology
iceRheol = 2; %iceRheol: ice rheology choice []
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

%%%%%     Pseudo-Initial Conditions     %%%%%
u(:,1) = zeros(m*n,1); %u: velocity[m/s]
T(:,1) = 260.15*ones(m*nT,1); %T: temperature[K]
mu(:,1) = 1e14*ones(m*n,1); %mu: effective viscosity[Pa s]
k(:,1) = k1*exp(-k2*1e-3.*T); %k: thermal conductivity[W/m*K]
c(:,1) = c1 + c2.*T; %c: specific heat[J/kg*K]

%%%%%     Basal Strength Evolution     %%%%%
drivingRamp = 1;
enforcedMargin = 17e3; %enforcedMargin: margin location [m]
tau_base = -9.00e3 + (kk-1)*.15e3 - Z*rho*g*sin(alpha)*(((y)./(60*Z)).^50) ...
                   - (kk-1)*2e3*gaussmf(y,[kk*0.1e3,30e3]);
% tau_base = -2.96e3 - Z*rho*g*sin(alpha)*(((y)./(5*Z)).^10) ...
%                    - (kk-1)*1e3*gaussmf(y,[0.25e3,3e3]);
%tau_base: basal yield stress [Pa]


%% %% %% %% %%     Solution     %% %% %% %% %%

%%%%%     Steady-State Initial Conditions     %%%%%
fprintf(fID,'*****     Time: %d     *****\n',t(1));
[u(:,1),T(:,1),mu(:,1),k(:,1),c(:,1),MTP] = ...
    marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,u(:,1),T(:,1),mu(:,1),...
                k(:,1),c(:,1),tau_base(:,1),enforcedMargin(:,1),north,...
                south,east,west,northT,southT,eastT,westT,iceRheol,...
                marginSolveType,rho,g,alpha,drivingRamp(1),k1,k2,c1,c2,...
                T_m,T_atm,G_base,w,v,omegaM,omegaT,tol,0);

%%%%%     Time Loop     %%%%%
for nt = 2:length(t)
    fprintf(fID,'*****     Time: %d     *****\n',t(nt));
    [u(:,nt),T(:,nt),mu(:,nt),k(:,nt),c(:,nt),MTP] = ...
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

% subplot(3,1,1)
% plotSurfaceVelocity(m,n,y,u);
% subplot(3,1,2)
% plotBasalStress(0,0,m,y,dz,Z,rho,g,alpha,tau_base,mu,u,0)
% subplot(3,1,3)
% plotTemperatureField(m,nT,y,zT,T,T_m)

%Compute integrated stresses
u_surf_diff2 = gradient(Z*mean(reshape(mu,m,n),2).*...
               gradient(mean(reshape(u,m,n),2))/dy)/dy;
for i = 1:m
    tau_base_int(i) = trapz(-tau_base(1:i))*dy;
    tau_grav_int(i) = trapz(rho*g*Z*sin(alpha)*ones(i,1))*dy;
    tau_lat_int(i) = trapz(-u_surf_diff2(1:i))*dy;
end
%Plot integrated stresses
hFig = figure(2);
set(hFig, 'Position', [0 0 1700 400])
subplot(1,3,1)
plot(MTP/1e3*ones(2,1),[0;1],'k--')
hold on
plot(y/1e3,tau_grav_int/(rho*g*Z*sin(alpha)*Y),'LineWidth',3)
plot(y/1e3,tau_base_int/(rho*g*Z*sin(alpha)*Y),'LineWidth',3)
plot(y/1e3,tau_lat_int/(rho*g*Z*sin(alpha)*Y),'LineWidth',3)
plot(y/1e3,(tau_lat_int+tau_base_int)/(rho*g*Z*sin(alpha)*Y),'LineWidth',3)
legend('Margin Position','\tau_{grav}','\tau_{base}','\tau_{lat}',...
       '\tau_{base}+\tau_{lat}','Location','NW')
axis([0 Y/1e3 0 1])
xlabel('Distance from stream center [km]')
ylabel('Normalized integrated stress')
%Plot basal stresses
subplot(1,3,2)
plotBasalStress(0,0,m,y,dz,Z,rho,g,alpha,tau_base,mu,u,MTP)
hold on
plot(MTP/1e3*ones(2,1),[0;100],'k--')
subplot(1,3,3)
plotSurfaceVelocity(m,n,y,u);
setFontSize(13)
anim(kk) = getframe(hFig);
close all

end