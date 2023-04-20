function[] = plotStressStrain(z,y,u,mu,m,n,dy,dz)
%Plot the stress and strain fields
[tau_E,epsilon_E] = effectiveStress(u,mu,m,n,dy,dz);

ax3 = subplot(311);
[C,Ch] = contourf(y/1000,z,reshape(u*3.154e7,m,n)');
set(Ch,'LineColor','none')
colormap(ax3, cmocean('Tempo'))
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')
title('Speed [m/yr]')


ax1 = subplot(312);
[C,Ch] = contourf(y/1000,z,reshape(log10(epsilon_E),m,n)');
colormap(ax1, cmocean('Speed'))
set(Ch,'LineColor','none')
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')
title('Log10 Effective Strain')

ax2 = subplot(313);
[C,Ch] = contourf(y/1000,z,reshape(log10(tau_E),m,n)');
set(Ch,'LineColor','none')
colormap(ax2, cmocean('Turbid'))
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')
title('Log10 Effective Stress')
end

%Template
% [C,Ch] = contourf(y/1000,zT,reshape(T-T_m,m,nT)');
% clabel(C,Ch);
% xlabel('Distance from Stream Center [km]')
% ylabel('Distance from Base [m]')