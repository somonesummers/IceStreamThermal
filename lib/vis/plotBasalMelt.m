function[] = plotBasalMelt(y,m,n,nT,dz,T,T_m,tau_E,epsilon_E,rho,L,G_base,u,tau_base)
%Plot the basal melt flux at the base

b_dot = meltwaterFlux(m,nT,dz,T,T_m,[zeros(m*(nT-n),1);tau_E],...
                      [zeros(m*(nT-n),1);epsilon_E],G_base,rho,L,u,tau_base);

plot(y/1e3,-b_dot*1000*3.15569e7)
xlabel('Distance from Stream Center [km]')
ylabel('Melt Rate [mm/yr]')

end

