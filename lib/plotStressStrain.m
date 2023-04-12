function[] = plotStressStrain(z,y,u,mu,m,n,dy,dz)
%Plot the stress and strain fields
[tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz);

ax1 = subplot(211);
[C,Ch] = contourf(y/1000,z,reshape(epsilon_E,m,n)');
colormap(ax1, cmocean('Speed'))
set(Ch,'LineColor','none')
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')
title('Effective Strain')

ax2 = subplot(212);
[C,Ch] = contourf(y/1000,z,reshape(tau_E,m,n)');
set(Ch,'LineColor','none')
colormap(ax2, cmocean('Turbid'))
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')
title('Effective Stress')
end

%Template
% [C,Ch] = contourf(y/1000,zT,reshape(T-T_m,m,nT)');
% clabel(C,Ch);
% xlabel('Distance from Stream Center [km]')
% ylabel('Distance from Base [m]')