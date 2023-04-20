function [mu] = iceRheology_goldsby2000(T,tau_E,m,n,rho,g,z)
%Compute non-Newtonian ice rheology from Goldsby

%%%%%     Physical Constants     %%%%%
P = reshape(rho*g*repmat(flipud(z)',m,1),m*n,1); %P: hydrostatic pressure[Pa]
p0 = 7e-8; %p0: pressure heating coef[K/Pa]
R = 8.314; %R: gas constant[J/K*mol]
%Diffusional Rheology
Omega = 3.27e-29; %Omega: molecular volume[m^3]
kB = 1.38e-23; %k_B: Boltzmann constant[m^2*kg/s*K]
d = 5e-3; %d: grain size[m]
B = 9.1e-4; %B: exponential prefactor[m^2/s]
Q_diff = 59.4e3; %Q1: diffusional activation energy[J/mol]
%Basal Rheology
A_basal = 2.19e-7; %preexponential constant[Pa^2.4]
Q_basal = 60e3;
%GBS Rheology
A1_gbs = 4.0e-23; %A_gbs: preexponential constant[m^1.4/s*Pa^1.8]
A2_gbs = 4.0e-23; %A_gbs: preexponential constant[m^1.4/s*Pa^1.8]
% A1_gbs = 6.18e-14; %A_gbs: preexponential constant[m^1.4/s*Pa^1.8]
% A2_gbs = 8.24e16; %A_gbs: preexponential constant[m^1.4/s*Pa^1.8]
Q1_gbs = 49e3; %Q1_gbs: lower activation energy[J/mol]
Q2_gbs = 197e3; %Q2_gbs: higher activation energy[J/mol]
Tstar_gbs = 255; %Tstar_gbs: activation threshold[K] (-18C)
Q_gbs_h = zeros(m*n,1); %Q_gbs_h: melting point adjusted activation energy[J/mol]
A_gbs  = zeros(m*n,1); %Q_gbs_h: melting point adjusted  preexponential constant [m^1.4/s*Pa^1.8]
%Dislocation Rheology
% A_disl = 3.162e-30; %A_disl: preexponential constant[1/s*Pa^4] Peltier 2000
% A_disl = 3.162e-30; %A_disl: preexponential constant[1/s*Pa^4] (Qi 2021)
A_disl = 6.0e-33; %A_disl: preexponential constant[1/s*Pa^4] Paul's tuning
Q1_disl = 64e3; %Q1_disl: lower activation energy[J/mol]
Q2_disl = 220e3; %Q2_disl: higher activation energy[J/mol]
Tstar_disl = 255; %Tstar_gbs: activation threshold[K] (-15C)
Q_disl_h = zeros(m*n,1); %Q_disl_h: melting point adjusted activation energy[J/mol]

%Compute spatially varying activation energies
T_h = T + p0*P; %T_h: melting point adjusted temperature[K]
Tstar_gbs_h = Tstar_gbs*ones(m*n,1) + p0*P;
Tstar_disl_h = Tstar_disl*ones(m*n,1) + p0*P; %Tstar_h: adjusted activation threshold[K]
A_gbs(heaviside(T_h-Tstar_gbs_h)==0) = A1_gbs;
A_gbs(heaviside(T_h-Tstar_gbs_h)>0.5) = A2_gbs;
Q_gbs_h(heaviside(T_h-Tstar_gbs_h)==0) = Q1_gbs;
Q_gbs_h(heaviside(T_h-Tstar_gbs_h)>=0.5) = Q2_gbs;
Q_disl_h(heaviside(T_h-Tstar_disl_h)==0) = Q1_disl;
Q_disl_h(heaviside(T_h-Tstar_disl_h)>=0.5) = Q2_disl;

eps_diff = ((42*Omega*B)./(kB*(d^2)*T_h)).*exp(-Q_diff./(R*T_h));%*tau_ij
eps_basal = A_basal.*exp(-Q_basal./(R*T_h)).*tau_E.^1.4;%*tau_ij
eps_gbs = (A_gbs./(d^1.4)).*exp((-Q_gbs_h./R).*((1./T_h)-(1./Tstar_gbs_h))).*tau_E.^0.8;
eps_disl = A_disl.*exp((-Q_disl_h./R).*((1./T_h)-(1./Tstar_disl_h))).*tau_E.^3;

% mu = 0.5*(eps_diff + eps_gbs + eps_disl).^-1;
mu = 0.5*(eps_diff + (eps_basal.*eps_gbs)./(eps_basal + eps_gbs) + eps_disl).^-1;

end

%Visualization of each contribution
% subplot(1,4,1)
% contourf(y,z,reshape(log10(eps_diff),m,n)')
% colorbar
% title('Diffusional Creep')
% subplot(1,4,2)
% contourf(y,z,reshape(log10(eps_basal),m,n)')
% colorbar
% title('Basal Slip')
% subplot(1,4,3)
% contourf(y,z,reshape(log10(eps_gbs),m,n)')
% colorbar
% title('GBS Slip')
% subplot(1,4,4)
% contourf(y,z,reshape(log10(eps_disl),m,n)')
% colorbar
% title('Dislocation Creep')