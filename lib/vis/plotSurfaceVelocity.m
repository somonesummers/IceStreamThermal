function[] = plotSurfaceVelocity(m,n,y,u)
%Plot surface velocity and compare to Echelmeyer and Harrison (1999)

plot(y(1:m)/1000,3.1557e7*u(m*n-m+1:m*n),'k')
ylabel('Surface Velocity [m/yr]')
xlabel('Distance from Stream Center [km]')

end

%Original template
% load('data/whillansSurfVel.mat')
% plot(y(1:m)/1000,3.1557e7*u(m*n-m+1:m*n),y_ob+obOffset,u_ob,'-o')
% axis([0,22,0,500])