function [ LOSinfo, time, z ] = calc_LoS_OSIRIS( tag )
%CALC_LOS_GBS Calculates 'line of sight' information OSIRIS - use only
%tangent point!!
%   INPUT:
%       tag ('NO2' or 'O3'): select tracegas file
%
%   OUTPUT:
%       LOSinfo.Lat: latitude of point with altitude z along LoS (length(z) x length(time))
%       LOSinfo.Lon: longitude of point with altitude z along LoS (length(z) x length(time))
%       LOSinfo.z, z: altitude levels
%       time: structure with date and time info

%% 

switch tag
    case 'NO2'
        load('/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/OSIRIS_v6p00_NO2_table_12-32km.mat');
    case 'O3'
        load('/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/OSIRIS_v5p10_O3_table_14-52km.mat');
    otherwise
        error('First input can only be NO2 or O3')
end

z=fliplr([0.795, [1.5:1:59.5]]);

% fil LoS structure
LOSinfo=struct;
LOSinfo.z=z;

LOSinfo.Lat=repmat(osiris.lat',length(z),1);
LOSinfo.Lon=repmat(osiris.lon',length(z),1);

% create time structure
time_array=mjd2k_to_date(osiris.mjd2k);

time=struct;

time.year=year(time_array);
time.month=month(time_array);
time.day=day(time_array);
time.hour=hour(time_array);
time.min=minute(time_array);
time.sec=second(time_array);

time.UTC=zeros(size(time.year));

end


