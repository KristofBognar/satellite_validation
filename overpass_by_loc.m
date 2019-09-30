function [ ace, osiris ] = overpass_by_loc( loc_in )
%OVERPASS_BY_LOC find satellite overpasses for given locations
%   location can be passed either as a structure with lat, lon and alt
%   fields, or as a string of the location name

%% location data


if isstruct (loc_in)
    loc=loc_in;
else
    switch loc_in
        case 'alert'
            loc.lat=82.5;
            loc.lon=-62.35;
            loc.alt=10;
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
    end
end
%% ACE data

data_dir = '/home/kristof/work/satellite_validation/ACE-FTS_data/';

% use ozone data
load([data_dir 'ACE_v3p6_O3.mat']);

ace=table;

% filter by geolocation
% only include measurements within 500 km of location
range_km=dist_to_loc(tanstruct.lat_tangent,tanstruct.lon_tangent,loc);

% filter by distance
goodind=find(range_km<=500);

ace.lat=tanstruct.lat_tangent(goodind)';
ace.lon=tanstruct.lon_tangent(goodind)';


%% OSIRIS data

cur_dir=(pwd);
data_dir='/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/';

cd([data_dir 'O3']);

tmp=dir('*.nc');
flist={tmp.name};
flist(1:2)=[];

osiris=[];   

% loop through files and read in coordinates
for i=1:length(flist)
    
    % load coordinates
    os_lat=ncread(flist{i},'latitude');
    os_lon=ncread(flist{i},'longitude');

    % filter by geolocation (only keep measurements within 500km of PEARL)
    range_km=dist_to_loc(os_lat,os_lon,loc);
    goodind=find(range_km<=500);

    if ~isempty(goodind)
        osiris=[osiris; [double(os_lat(goodind)), double(os_lon(goodind))]];
    end

end

osiris= array2table(osiris,'VariableNames',{'lat','lon'});

cd(cur_dir)

end

