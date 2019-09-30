% plot measurement locations on a map

% gridm: Toggle and control display of graticule lines (lat lon grid)
% setm: Set properties of map axes and graphics objects
% mlabel: Toggle and control display of meridian labels

% axes properties:
% https://www.mathworks.com/help/map/ref/mapaxes-properties.html


eureka = [80.05, -86.42];
load coast


figure

hold on

ax = worldmap([70,85], [-120,-50]);
geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])

% plotm(ace_fts.lat,ace_fts.lon,'k.','markersize',8)
% plotm(osiris.lat,osiris.lon,'k.','markersize',8)
% osiris.lat=double(osiris.lat);
% osiris.lon=double(osiris.lon);

% % location.latitude=80.053;
% % location.longitude=-86.416;
% % location.altitude=610;
% % 
% % load('/home/kristof/work/GBS/VCD_results/UT-GBS_O3_VCD_all.mat')
% % 
% % LOSinfo = calc_LoS_GBS_interp( 'O3_VIS', location, reanalysis.sza, reanalysis.saa+180;, z);
% %
% % plotm(LOSinfo.Lat(31,:),LOSinfo.Lon(31,:),'k.','markersize',8)


plotm(80.05, -86.42, 'rp','markerfacecolor','r','markersize',10)


gridm('MLineLocation',10)
mlabel('MLabelLocation',10)

