function [ time, location, sza_out, saa_out, z ] = sza_saa_for_LoS_meas_time( instr_name, yr_out )
%SZA_SAA_FOR_LOS Summary of this function goes here
%   Detailed explanation goes here

%% setup

do_saoz=true;

% location info
location.latitude=80.053;
location.longitude=-86.416;
location.altitude=610;

z=fliplr([0.795, [1.5:1:59.5]]);

% data folder
data_dir='/home/kristof/work/DMP/all_GBS_meas_times/';

% file naming
if strcmp(instr_name,'O3_VIS')        
    str_tmp='O3';        
elseif strcmp(instr_name,'NO2_VIS')        
    str_tmp='NO2';        
elseif strcmp(instr_name,'NO2_UV')        
    str_tmp='NO2_UV';        
    do_saoz=false;
else error('Valid trace gas inputs are O3_VIS, NO2_VIS, and NO2_UV')
end

% check if only subset of years is needed
if length(yr_out)==5 && yr_out(5)=='+' % input year and later
    
    do_loop=1;
    yr_end=year(datetime(now,'convertfrom','datenum'));
    yr_start=str2double(yr_out(1:4));
    
elseif length(yr_out)==9 && yr_out(5)=='-' % year range
    
    do_loop=1;
    yr_end=str2double(yr_out(6:9));
    yr_start=str2double(yr_out(1:4));
    
else % one year or all years
    do_loop=0;
end

%% merge UT-GBS, PEARL-GBS, and SAOZ datasets for each species
meas_all=[];

if ~do_loop % one year or all years: single file per instrument
    
    %load GBS data and save time and solar position (use table fields,
    %location of columns is inconsistent)
    try 
        load([data_dir 'UT-GBS_' str_tmp '_VCD_' yr_out '_unfiltered.mat']);
        meas_all=[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]; 
    end

    try 
        load([data_dir 'PEARL-GBS_' str_tmp '_VCD_' yr_out '_unfiltered.mat']);
        meas_all=[meas_all;[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]];
    end

else % multiple years
        
    for i=yr_start:yr_end
        
        try 
            load([data_dir 'UT-GBS_' str_tmp '_VCD_' num2str(i) '_unfiltered.mat'])
            meas_all=[meas_all;[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]]; 
        end

        try 
            load([data_dir 'PEARL-GBS_' str_tmp '_VCD_' num2str(i) '_unfiltered.mat'])
            meas_all=[meas_all;[reanalysis.sza,reanalysis.saa,reanalysis.mjd2k]];
        end
        
    end
    
end

% convert to N=0 for saa (SolarAzEl returns saa with N=0)
meas_all(:,2)=meas_all(:,2)+180;

if do_saoz

    % load SAOZ data
    load([data_dir 'saoz_' lower(str_tmp) '.mat'])

    if strcmp(yr_out,'all')
        ind=true(size(saoz.year));
    else
        if length(yr_out)==5 && yr_out(5)=='+' % input year and later
            ind=(saoz.year>=str2double(yr_out(1:4)));
        elseif length(yr_out)==4 % input year only
            ind=(saoz.year==str2double(yr_out));
        end
    end
    
    % add to merged table
    meas_all=[meas_all;[saoz.sza(ind),saoz.saa(ind),saoz.mjd2k(ind)]];
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

