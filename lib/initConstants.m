function [rho,g,alpha,L,G_base,k1,k2,c1,c2,T_m,T_atm,...
          dy,dz,y,z,zT,nT,t,north,south,east,west,...
          northT,southT,eastT,westT,w,v] = initConstants(m,n,Y,Z,Z_till,endT,dt,v_0,a)
%Initialize constants

%%%%%     Physical Constants     %%%%%
%Mechanical
rho = 917; %rho: density[kg/m^3]
% rho_w = 1000; %rho_w: density of water[kg/m^3]
% mu_w = 1.787e-3; %mu_w: dynamic viscosity of water[Pa*s]
g = 9.81; %g: gravity[m/s^2]
alpha = 0.0012; %alpha: inclination angle[rad]
%Thermal
L = 335e3; %L: latent heat per unit mass[J/kg]
G_base = 80e-3; %46.5e-3; %G_base: geothermal heat flux[W/m^2]
k1 = 9.828; %k1: conductivity preexponential[W/m*K]
k2 = 5.7; %k2: conductivity postexponential[1/K]
c1 = 152.5; %c1: specific heat forefactor[J/kg*K]
c2 = 7.122; %c2: specific heat forefactor[J/kg*K^2]
T_m = 273.15; %T_m: melting temperature[K]
T_atm = 247.15; %T_atm: atmospheric temperature[K]

%%%%%     Preallocate     %%%%%
dy = Y/(m-1);
dz = Z/(n-1); %dy,dz: triangulation width[m]
y = (0:dy:Y)';
z = (0:dz:Z)'; %y,z: triangulation axes
zT = (-Z_till:dz:Z)'; %z_till: thermal triangulation with till
nT = length(zT); %nT: vertical dimension of thermal domain
t = (0:dt:endT)'; %t: time discretization
north = (m*(n-1)+1:m*n)';
south = (1:m)';
east = (m:m:n*m)';
west = (1:m:m*(n-1)+1)'; %north,...: boundary node numbers on each face
northT = (m*(nT-1)+1:m*nT)';
southT = (1:m)';
eastT = (m:m:nT*m)';
westT = (1:m:m*(nT-1)+1)'; %north,...: thermal boundary node numbers on each face
v = v_0*y; %v: horizontal advection [m/s]
w = (-a*3.169e-8/Z)*z; %w: vertical advection [m/s]

end

