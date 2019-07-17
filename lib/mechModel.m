function [u,tau_E,epsilon_E] = mechModel(m,n,y,dy,dz,Z,mu,tau_base,H_m,...
                                         north,south,east,west,rho,g,...
                                         alpha,drivingRamp)
%Solves mechanical model

%%%%%     Build System     %%%%%
M = secondOrderLaplacian2D(m,n,dy,dz,mu,ones(m*n,1));
b = -rho*g*sin(alpha)*ones(m*n,1)*drivingRamp; %b: body force[kg/m s^2]

%%%%%     Apply BCs     %%%%%
%North: stress-free, South: Robin, East,West: symmetry
[M,b] = applyMechBC(M,b,mu,H_m,tau_base,north,south,east,west,m,dy,dz);

%%%%%     Solve system     %%%%%
u = M\b;

%%%%%     Compute effective stress     %%%%%
[tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz); %tau_E: effective stress[Pa]

end

