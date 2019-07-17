function [mu,resM] = mechRelax(m,n,nT,z,T,mu,tau_E,rho,g,omegaM,iceRheol)
%Perform under-relaxation on the viscosity

mu_old = mu;
if any(iceRheol==[2,3])
    mu = iceRheology_suckale(T(m*(nT-n)+1:end),tau_E,m,n,rho,g,z);
elseif iceRheol==4
    mu = iceRheology_goldsby(T(m*(nT-n)+1:end),tau_E,m,n,rho,g,z);
end
mu = omegaM*mu + (1-omegaM)*mu_old;
resM = norm((mu_old-mu)./mu_old);

end

