function[] = plotBuriedCrevasses(m,n,dy,y,Z,u,t,a)
%Plot buried crevasses advected downwards by accumulation

w =@(z) (-a*3.169e-8/Z)*z;
for i = 1:length(t)
    indL = find(gradient(3.1557e7*u(m*n-m+1:m*n,i))/dy <= -0.02,1,'first');
    indR = find(gradient(3.1557e7*u(m*n-m+1:m*n,i))/dy <= -0.02,1,'last');
    pos = [y(indL);y(indR)];
    depth = Z;
    for ii = 1:1000
        depth = depth + w(depth)*(t(end)-t(i))/1000;
    end
    plot(pos/1e3,depth*(pos./pos),'r','LineWidth',2)
    hold on
end

end