function [ range_km ] = dist_to_loc( lat_arr, lon_arr, location )
%[ range_km ] = dist_to_loc( lat_arr, lon_arr )
%   Calculates distance of Coordinates to any location

% create matching arrays
loc_lat=ones(size(lat_arr))*location.lat;
loc_lon=ones(size(lon_arr))*location.lon;

% calculate arc lengths
[arclen,~]=distance(loc_lat,loc_lon, lat_arr,lon_arr);

% calculate distances
R_e = 6378.1;
range_km = arclen * (pi/180) * R_e;


end

