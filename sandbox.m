% Site 1 -77.3281  -100.0582
% Depth 3 km
% Slope = 0.0012
% speed @ 25km in-stream 50 m/yr
xi = 4000;
phi = pi/3;
[x,y] = ll2ps(-77.3281,-100.0582);
pathx = [x - xi*cos(phi), x + xi*cos(phi)];
pathy = [y - xi*sin(phi), y + xi*sin(phi)];
dist = pathdistps(pathx,pathy,'m');

hice = bedmap2_interp(pathx,pathy,'surface');
slope = (hice(2)-hice(1))/dist(2);



% Site 2 -76.4085  -103.4856
% Depth 1.8 km
% Slope 0.0037
% speed @ 25km in-stream 100 m/yr
xi = 4000;
phi2 = pi/6;
[x,y] = ll2ps(-76.4085,-103.4856);
pathx2 = [x - xi*cos(phi2), x + xi*cos(phi2)];
pathy2 = [y - xi*sin(phi2), y + xi*sin(phi2)];
dist2 = pathdistps(pathx2,pathy2,'m');

hice2 = bedmap2_interp(pathx2,pathy2,'surface');
slope2 = (hice2(2)-hice2(1))/dist2(2);


%plot both
figure
measures('speed','thwaites glacier','scalelim',[1 600],'mapwidth',800);
[lat,lon] = ps2ll(pathx,pathy);
[lat2,lon2] = ps2ll(pathx2,pathy2);
plotm(lat,lon,'k','LineWidth',5)
plotm(lat2,lon2,'k','LineWidth',5)
plotm(-77.3281,-100.0582,'rp')
plotm(-76.4085,-103.4856,'rp')