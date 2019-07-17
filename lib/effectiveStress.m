function [tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz)
%Compute effective stress

[u_y,u_z] = gradient(reshape(u,m,n)',dy,dz);
u_y = reshape(u_y',m*n,1);
u_z = reshape(u_z',m*n,1);
tau_E = sqrt((mu.*u_y).^2 + (mu.*u_z).^2); %tau_E: effective stress[Pa]
epsilon_E = sqrt((u_y/2).^2 + (u_z/2).^2); %epsilon_E: effective strain rate[1/s]

end

