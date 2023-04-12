function[] = plotMuField(m,n,y,z,mu)
%Plot the temperature field

colormap parula
[C,Ch] = contourf(y/1000,z,log10(reshape(mu,m,n)'));
set(Ch,'LineColor','none')
colorbar
xlabel('Distance from Stream Center [km]')
ylabel('Distance from Base [m]')

end

%Template
% [C,Ch] = contourf(y/1000,zT,reshape(T-T_m,m,nT)');
% clabel(C,Ch);
% xlabel('Distance from Stream Center [km]')
% ylabel('Distance from Base [m]')