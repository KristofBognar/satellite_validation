function [ range_km ] = dist_to_PEARL( lat_arr, lon_arr )
%[ range_km ] = dist_to_PEARL( lat_arr, lon_arr )
%   Calculates distance of Coordinates to PEARL Ridge Lab


% location info for PEARL
pearl.lat=80.053;
pearl.lon=-86.416;
pearl.alt=610;

% create matching arrays
pearl_lat=ones(size(lat_arr))*pearl.lat;
pearl_lon=ones(size(lon_arr))*pearl.lon;

% calculate arc lengths
[arclen,~]=distance(pearl_lat,pearl_lon, lat_arr,lon_arr);

% calculate distances
R_e = 6378.1;
range_km = arclen * (pi/180) * R_e;


end

