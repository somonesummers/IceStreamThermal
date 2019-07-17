function [b_dot] = meltwaterFlux(m,nT,dz,T,T_m,tau_E,epsilon_E,...
                                  G_base,rho,L,u,tau_base)
%Construct meltwater flux at the bed

%Englacial melt(shear heating)
H = heaviside(T-T_m);
H(H==0.5) = 1;
S = H.*(2*tau_E.*epsilon_E);
b_e = dz*(S(1:m)./2 + sum(reshape(S(m+1:end-m),m,nT-2),2) + S(end-m+1:end)./2);

%Basal abalation
H_BC = heaviside(T(m+1:2*m)-T_m);
H_BC(H_BC == 0.5) = 1;
H_BC(1:find(H_BC==1,1)) = 1;
b_b = H_BC.*(G_base - u(1:m).*tau_base);

%Basal meltwater flux
b_dot = -(b_e + b_b)/(rho*L);

end

