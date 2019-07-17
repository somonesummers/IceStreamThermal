clear all; clc; path(pathdef);
addpath('~/Documents/Research/iceStreamResearch/mapData/modismoa_v4/')
addpath('~/Documents/Research/iceStreamResearch/mapData/bedmap2_tiff/')

%Plot sample lines on MODIS image
subplot(2,3,[1,2,4,5])
modismoa(-82.6,-146.7,350,'uhc')
plotm([-82.0;-83.3],[-150.8;-150.0],'k','lineWidth',5)
plotm([-82.8;-82.5],[-150.8;-145.0],'k','lineWidth',5)

%Plot across-stream thickness
[lat,lon,z] = bedmap2_data('thickness','resolution','2 km');
lat = reshape(lat,length(lat)^2,1);
lon = reshape(lon,length(lon)^2,1);
z = reshape(z,length(z)^2,1);
F = TriScatteredInterp(lat,lon,z,'natural');
line_lat = linspace(-82.0,-83.3,100)';
line_long = linspace(-150.8,-150.0,100)';
line_dist = pathdist(line_lat,line_long);
line_thickness = F(line_lat,line_long);
subplot(2,3,3)
plot(line_dist/1e3,line_thickness)

%Plot along-stream surface/bed elevation
clear all
[lat,lon,z] = bedmap2_data('bed','resolution','2 km');
lat = reshape(lat,length(lat)^2,1);
lon = reshape(lon,length(lon)^2,1);
z = reshape(z,length(z)^2,1);
F = TriScatteredInterp(lat,lon,z,'natural');
line_lat = linspace(-82.5,-82.8,100)';
line_long = linspace(-145.0,-150.8,100)';
line_dist = pathdist(line_lat,line_long);
line_surface = F(line_lat,line_long);
subplot(2,3,6)
plot(line_dist/1e3,line_surface,line_dist/1e3,-0.6141*(line_dist/1e3)+205)