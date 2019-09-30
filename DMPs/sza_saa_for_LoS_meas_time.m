function [ time, location, sza_out, saa_out, z ] = sza_saa_for_LoS_meas_time( instr_name )
%SZA_SAA_FOR_LOS Summary of this function goes here
%   Detailed explanation goes here

%% setup

% location info
location.latitude=80.053;
location.longitude=-86.416;
location.altitude=610;

z=fliplr([0.795, [1.5:1:59.5]]);

%% merge UT-GBS, PEARL-GBS, and SAOZ datasets for each species

do_saoz=true;

if strcmp(instr_name,'O3_VIS')        
    str_tmp='O3';        
elseif strcmp(instr_name,'NO2_VIS')        
    str_tmp='NO2';        
elseif strcmp(instr_name,'NO2_UV')        
    str_tmp='NO2_UV';        
    do_saoz=false;
else error('Valid trace gas inputs are O3_VIS, NO2_VIS, and NO2_UV')
end

%load GBS data and save time and solar position (use table fields,
%location of columns is inconsistent)
load(['/home/kristof/work/DMP/all_GBS_meas_times/UT-GBS_' str_tmp '_VCD_all_unfiltered.mat'])
meas_all=[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]; 

load(['/home/kristof/work/DMP/all_GBS_meas_times/PEARL-GBS_' str_tmp '_VCD_all_unfiltered.mat'])
meas_all=[meas_all;[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]];

% convert to N=0 for saa (SolarAzEl returns saa with N=0)
meas_all(:,2)=meas_all(:,2)+180;

if do_saoz
    % load SAOZ data
    load(['/home/kristof/work/SAOZ/saoz_' lower(str_tmp) '.mat'])

    % add to merged table
    meas_all=[meas_all;[saoz.sza,saoz.saa,saoz.mjd2k]];
end


% sort by time
meas_all=sortrows(meas_all,3);

% convert back to table
meas_all=array2table(meas_all,'variablenames',{'sza','saa','mjd2k'});


%% output
% get date
time_array=mjd2k_to_date(meas_all.mjd2k);

time=struct;

time.year=year(time_array);
time.month=month(time_array);
time.day=day(time_array);
time.hour=hour(time_array);
time.min=minute(time_array);
time.sec=second(time_array);

time.UTC=zeros(size(time.year));

% get sza, saa
sza_out=meas_all.sza;
saa_out=meas_all.saa;

end

