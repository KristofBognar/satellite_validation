load coast

% figure()
% 
% set(gcf, 'Position', [100, 100, 1000, 600]);
% 
% ax = worldmap('world');
% geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])
% 
% plotm(tanstruct.lat_tangent(1:5000),tanstruct.lon_tangent(1:5000),...
%       '.','color',[0,0,0]+0.2,'markersize',8)
% plotm(80.053, -86.416, 'yp','markerfacecolor','y','markersize',13)
% 
% mlabel('off')
% plabel('off')
% 
% 
figure()

set(gcf, 'Position', [100, 100, 600, 600]);

ax = worldmap([65,90], [-180,180]);
geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])

plotm(tanstruct.lat_tangent(1:10000),tanstruct.lon_tangent(1:10000),...
      '.','color',[0,0,0]+0.4,'markersize',9)
plotm(80.053, -86.416, 'rp','markerfacecolor','r','markersize',13)

mlabel('off')
plabel('off')

% figure()

R_e = 6378.1;
arclen_500 = (500/R_e) * (180/pi);
arclen_1000 = (1000/R_e) * (180/pi);

i=1;
for az=1:0.1:360
    [lat_500(i),lon_500(i)]=reckon(80.053,-86.416,arclen_500,az);
    [lat_1000(i),lon_1000(i)]=reckon(80.053,-86.416,arclen_1000,az);
    i=i+1;
end

% ax = worldmap([65,85], [-136.4,-36.4]);
% geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])

% plotm(80.053,-86.416, 'rp','markerfacecolor','r','markersize',12)
plotm(lat_500,lon_500,'r','linewidth',2)
plotm(lat_1000,lon_1000,'r','linewidth',2)

textm(lat_500(2700),lon_500(2700)-5,'500 km','fontsize',14,'color','r',...
    'fontweight','bold','Horizontalalignment','center')
textm(lat_1000(2700),lon_1000(2700)-5,'1000 km','fontsize',14,'color','r',...
    'fontweight','bold','Horizontalalignment','center')




