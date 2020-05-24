function [time,location,sza,saa,z] = prepDMPinput_bruker_dat(yr_out)
%PREPDMPINPUT_BRUKER_DAT Summary of this function goes here
%   Detailed explanation goes here

%% output variables
time=struct;

location.latitude=80.053;
location.longitude=-86.416;
location.altitude=610;


%% load file with bruker measurement time info
load('/home/kristof/work/bruker/spDB_eur_2006_2020.mat')

z=fliplr(ALTITUDE);

% convert date and time to datetime (UTC time)
try
    tmp = strcat(data.Date, {' '}, data.Time);
    time_array = datetime(tmp,'InputFormat','yyyyMMdd HH:mm:ss');
catch
    error('Save spDB file as table, with date and time files as TEXT')
end

% check if only subset of years is needed
if strcmp(yr_out,'all')
    ind=true(size(time_array));
else
    if length(yr_out)==5 && yr_out(5)=='+' % input year and later
        ind=(time_array.Year>=str2double(yr_out(1:4)));
    elseif length(yr_out)==4 % input year only
        ind=(time_array.Year==str2double(yr_out));
    elseif length(yr_out)==9 && yr_out(5)=='-' % year range
        ind=(time_array.Year>=str2double(yr_out(1:4)) & time_array.Year<=str2double(yr_out(6:9)));
    end
end

% time output
time.year=year(time_array(ind));
time.month=month(time_array(ind));
time.day=day(time_array(ind));
time.hour=hour(time_array(ind));
time.min=minute(time_array(ind));
time.sec=second(time_array(ind));

time.UTC=zeros(size(time.year));

% solar position output
sza=data.SZen(ind);
saa=data.SAzm(ind);

end

