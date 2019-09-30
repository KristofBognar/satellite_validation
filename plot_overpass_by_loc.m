function [  ] = plot_overpass_by_loc( )
%PLOT_OVERPASS_BY_LOC plot number of ace/OSIRIS measurements at various locations
%

loc_list={'alert','nord','eureka','ny-alesund','resolute','utqiagvik'};

for i=1:length(loc_list)

    loc_in=loc_list{i};
    switch loc_in
        case 'alert'
            loc.lat=82.5;
            loc.lon=-62.35;
            loc.alt=10;
        case 'eureka'
            loc.lat=80.053;
            loc.lon=-86.416;
            loc.alt=610;
        case 'resolute'
            loc.lat=74.8;
            loc.lon=-94.83;
            loc.alt=10;
        case 'nord'
            loc.lat=81.72;
            loc.lon=-17.8;
            loc.alt=10;
        case 'ny-alesund'
            loc.lat=78.93;
            loc.lon=11.92;
            loc.alt=10;
        case 'cambridge'
            loc.lat=69.12;
            loc.lon=-105.05;
            loc.alt=10;
        case 'utqiagvik'
            loc.lat=71.29;
            loc.lon=-156.78;
            loc.alt=10;
        case 'thule'
            loc.lat=76.5;
            loc.lon=-68.8;
            loc.alt=10;
    end


    [ ace, osiris ] = overpass_by_loc( loc );

    figure(1)
    subplot(3,2,i)
    hold on
    map_tmp(loc,ace,loc_in)

    figure(2)
    subplot(3,2,i)
    hold on
    map_tmp(loc,osiris,loc_in)


end
end

function map_tmp(loc,sat, loc_in)
    
    load coast

    ax = worldmap([min(sat.lat)-1,max(sat.lat)+1], [min(sat.lon)-5,max(sat.lon)+5]);
    geoshow(ax, lat, long,'DisplayType', 'polygon', 'FaceColor', [0.7,0.7,0.7])

    plotm(sat.lat,sat.lon,'k.','markersize',8);

    plotm(loc.lat, loc.lon, 'rp','markerfacecolor','r','markersize',10)

    gridm('MLineLocation',20)
    mlabel('MLabelLocation',20)

    title({[upper(loc_in(1)), loc_in(2:end)] ,['n=' num2str(length(sat.lat))]})

end



