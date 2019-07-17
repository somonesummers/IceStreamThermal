function[] = plotTemperatureField(m,nT,y,zT,T,T_m)
%Plot the temperature field

load iceColorMap
colormap(iceColorMap)
[C,Ch] = contourf(y/1000,zT,reshape(T-T_m,m,nT)');
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