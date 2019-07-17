function [mu] = iceRheology(T,tau_E,m,n,rho,g,z)
%Compute non-Newtonian ice rheology from diffusional creep and Glen's Law

%%%%%     Physical Constants     %%%%%
P = reshape(rho*g*repmat(flipud(z)',m,1),m*n,1); %P: hydrostatic pressure[Pa]
%Diffusional Rheology
Omega = 3.27e-29; %Omega: molecular volume[m^3]
kB = 1.38e-23; %k_B: Boltzmann constant[m^2*kg/s*K]
d = 5e-3; %d: grain size[m]
B = 9.1e-4; %B: exponential prefactor[m^2/s]
Q1 = 59.4e3; %Q1: diffusional activation energy[J/mol]
R = 8.314; %R: gas constant[J/K*mol]
%Glen's Law Rheology
A = 3.5e-25; %A: preexponential constant[1/s*Pa^3]
E = 1; %E: enhancement factor[]
p0 = 7e-8; %p0: pressure heating coef[K/Pa]
Tstar = 263.15; %Tstar: activation threshold[K] (-10C)
Q2 = 60e3; %Q2: lower activation energy[J/mol]
Q3 = 115e3; %Q3: higher activation energy[J/mol]
Q_h = zeros(m*n,1); %Q_h: melting point adjusted activation energy[J/mol]

T_h = T + p0*P; %T_h: melting point adjusted temperature[K]
Tstar_h = Tstar*ones(m*n,1) + p0*P; %Tstar_h: adjusted activation threshold[K]
Q_h(heaviside(T_h-Tstar_h)==0) = Q2;
Q_h(heaviside(T_h-Tstar_h)==1) = Q3;

mu = 0.5*(((42*Omega*B)./(kB*(d^2)*T)).*exp(-Q1./(R*T))+...
          A*E*exp((-Q_h./R).*((1./T_h)-(1./Tstar_h))).*tau_E.^2).^-1;

end

