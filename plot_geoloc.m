function plot_geoloc()
%PLOT_GEOLOC plot geolocation of ACE, OSIRIS, GBS and Bruker measurements
%   Uses DMP data for ground-based instruments


% eureka coordinates
loc.lat=80.053;
loc.lon=-86.416;
loc.alt=610;

% load data
load('/home/kristof/work/satellite_validation/all_data_nosmooth.mat')
figure()

%% OSIRIS geolocations

subplot(221)
plot_tmp(osiris_o3,loc,'a) OS')

%% ACE geolocations

subplot(222)
plot_tmp(ace_fts_o3,loc,'b) AF')

%% Bruker geolocations

[~,~,~,bruker_pos.lat, bruker_pos.lon ]=...
    match_DMP_bruker(bruker_o3.fractional_time,bruker_o3.year,30);

subplot(223)

ind=find(bruker_o3.fractional_time<105 & bruker_o3.year<2015);

plot_tmp(bruker_pos,loc,'c) BK',ind)

% bruker_pos.lat=bruker_pos.lat(ind);
% bruker_pos.lon=bruker_pos.lon(ind);
% plot_tmp(bruker_pos,loc,'c) BK')

%% GBS geolocations

[~,~,~,gbs_pos.lat,gbs_pos.lon]=...
    match_DMP_DOAS( gbs_o3.fractional_time, gbs_o3.year, 'O3_VIS', 30 );

% ind=find(gbs_o3.day>55 & gbs_o3.day<105);
ind=find(gbs_o3.sza_min>85.5 & gbs_o3.sza_max-gbs_o3.sza_min>4);

% plot maps
subplot(224)
plot_tmp(gbs_pos,loc,'d) GV',ind)

% gbs_pos.lat=gbs_pos.lat(ind);
% gbs_pos.lon=gbs_pos.lon(ind);
% plot_tmp(gbs_pos,loc,'d) GV')

end

function plot_tmp(pos, loc, title_str, ind_in)

    load coast % coastline data for map
    
%     color_str='color',[0,0,0]+0.8

    ax = worldmap([70,85], [-126.4,-46.4]);
    geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])

    if nargin==4
        plotm(pos.lat,pos.lon,'.','color',[0,0,0]+0.3,'markersize',8);
        
        pos_tmp.lat=pos.lat(ind_in);
        pos_tmp.lon=pos.lon(ind_in);
        
        plotm(pos_tmp.lat,pos_tmp.lon,'k.','markersize',8);
        
    else
        plotm(pos.lat,pos.lon,'k.','markersize',8);
    end
    
    plotm(loc.lat, loc.lon, 'rp','markerfacecolor','r','markersize',12)
    
    gridm('MLineLocation',20)
    mlabel('MLabelLocation',20)

    title(title_str,'FontWeight','Normal')
%     tightmap on

end