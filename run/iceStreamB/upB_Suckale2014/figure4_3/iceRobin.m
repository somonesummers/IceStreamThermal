clear all; clc; path(pathdef);
addpath('../../../../lib/','../../../../lib/vis')
fID = fopen('log.txt','w');

%% %% %% %% %%     Initialization     %% %% %% %% %%

%%%%     Runtime variables     %%%%%
%Domain
m = 3401;
n = 101; %m,n: dimensions of triangulation
Y = 34e3;
Z = 1e3; %Y,Z: dimensions of domain[m]
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
marginSolveType = 2;
%1: self-consistent margin  %2: choose margin location
%Coupling
omegaT = 0.3; %omegaT: thermal relaxation parameter []
omegaM = 0.5; %omegaR: rheology relaxation parameter []
tol = 1e-3; %tol: rheological and thermal error tolerance []
%Initialize constants
[rho,g,alpha,G_base,k1,k2,c1,c2,T_m,T_atm,...
 dy,dz,y,z,zT,nT,t,north,south,east,west,...
 northT,southT,eastT,westT,w,v] = initConstants(m,n,Y,Z,Z_till,endT,dt,v_0,a);

%%%%%     Pseudo-Initial Conditions     %%%%%
u(:,1) = zeros(m*n,1); %u: velocity[m/s]
T(:,1) = 260.15*ones(m*nT,1); %T: temperature[K]
mu(:,1) = 1e14*ones(m*n,1); %mu: effective viscosity[Pa s]
k(:,1) = k1*exp(-k2*1e-3.*T); %k: thermal conductivity[W/m*K]
c(:,1) = c1 + c2.*T; %c: specific heat[J/kg*K]

%%%%%     Basal Strength Evolution     %%%%%
enforcedMargin = 17e3; %enforcedMargin: margin location [m]
tau_base = -[4.00e3*ones(floor(enforcedMargin*m/Y),1);...
             1e7*ones(m-floor(enforcedMargin*m/Y),1)];
%tau_base: basal yield stress [Pa]


%% %% %% %% %%     Solution     %% %% %% %% %%

%%%%%     Steady-State Initial Conditions     %%%%%
fprintf(fID,'*****     Time: %d     *****\n',t(1));
[u(:,1),T(:,1),mu(:,1),k(:,1),c(:,1)] = ...
    marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,u(:,1),T(:,1),mu(:,1),...
                k(:,1),c(:,1),tau_base(:,1),enforcedMargin(:,1),north,...
                south,east,west,northT,southT,eastT,westT,iceRheol,...
                marginSolveType,rho,g,alpha,k1,k2,c1,c2,T_m,T_atm,...
                G_base,w,v,omegaM,omegaT,tol,0);

%%%%%     Time Loop     %%%%%         
for nt = 2:length(t)
    fprintf(fID,'*****     Time: %d     *****\n',t(nt));
    [u(:,nt),T(:,nt),mu(:,nt),k(:,nt),c(:,nt)] = ...
        marginSolve(fID,m,n,nT,dy,y,dz,z,zT,Y,Z,dt,u(:,nt-1),T(:,nt-1),...
                    mu(:,nt-1),k(:,nt-1),c(:,nt-1),tau_base(:,nt),...
                    enforcedMargin(:,nt),north,south,east,west,northT,...
                    southT,eastT,westT,iceRheol,marginSolveType,rho,g,...
                    alpha,k1,k2,c1,c2,T_m,T_atm,G_base,w,v,omegaM,omegaT,tol,1);
    plotTemperatureField(m,nT,y,zT,T(:,end),T_m);
    getframe;
end


%% %% %% %% %%     Visualization     %% %% %% %% %%
save data.mat

% plotSurfaceVelocity(m,n,y,u);

% plotBasalStress(offset,eps,m,y,dz,Z,rho,g,alpha,tau_base,mu,u,marginPosition)

% plotMarginTemp(m,nT,y,zT,T,T_m,marginPosition)
