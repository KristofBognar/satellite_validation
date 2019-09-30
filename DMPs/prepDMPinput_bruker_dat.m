function [time,location,sza,saa,z] = prepDMPinput_bruker_dat()
%PREPDMPINPUT_BRUKER_DAT Summary of this function goes here
%   Detailed explanation goes here

%% output variables
time=struct;

location.latitude=80.053;
location.longitude=-86.416;
location.altitude=610;


%% load file with bruker measurement time info
load('/home/kristof/work/bruker/spDB_eur_2006_2017.mat')

z=fliplr(ALTITUDE);

% convert date and time to datetime (UTC time)
tmp = strcat(data.Date, {' '}, data.Time);

time_array = datetime(tmp,'InputFormat','yyyyMMdd HH:mm:ss');

% time output
time.year=year(time_array);
time.month=month(time_array);
time.day=day(time_array);
time.hour=hour(time_array);
time.min=minute(time_array);
time.sec=second(time_array);

time.UTC=zeros(size(time.year));

% solar position output
sza=data.sza;
saa=data.saa;

end

