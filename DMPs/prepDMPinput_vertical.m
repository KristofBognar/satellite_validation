function [time, LOSinfo] = prepDMPinput_vertical( yr_out )
%CALC_LOS_VERTICAL Vertical line-of-sight for DMPS: just returns the same
%coordinates with the required time and altitude dimensions

%%% set time in here OR load from file
% if set here, might just be hourly for the winter/spring

%%% select location in here OR load from file
% Ridge lab
loc_rl.latitude=80.053;
loc_rl.longitude=-86.416;
loc_rl.altitude=610;

% EWS
loc_ews.latitude=79.989;
loc_ews.longitude=-85.937;
loc_ews.altitude=10;

%%% set alt grid in here OR load from file

% select times
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

% get time structure
time=struct;

time.year=year(time_array(ind));
time.month=month(time_array(ind));
time.day=day(time_array(ind));
time.hour=hour(time_array(ind));
time.min=minute(time_array(ind));
time.sec=second(time_array(ind));

time.UTC=zeros(size(time.year));


% initialize output structure
LOSinfo=struct;
LOSinfo.z=z;
 
% assign fixed coordinates
% required shape is len(z) x len(time)
LOSinfo.Lat=repmat(location.latitude,length(z),length(time));
LOSinfo.Lon=repmat(location.longitude,length(z),length(time));



end

