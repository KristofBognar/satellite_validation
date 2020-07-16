function [time, LOSinfo, z] = prepDMPinput_vertical( loc_in, yr_out )
%CALC_LOS_VERTICAL Vertical line-of-sight for DMPS: just returns the same
%coordinates with the required time and altitude dimensions

%%% set alt grid 
z=fliplr([0:0.15:81]);

%%% select times
if strcmp(yr_out,'all')
    error('select start year')
    %yr1=;
    yr2=datetime(now,'convertfrom','datenum').Year;
else
    if length(yr_out)==5 && yr_out(5)=='+' % input year and later
        yr1=str2double(yr_out(1:4));
        yr2=datetime(now,'convertfrom','datenum').Year;
    elseif length(yr_out)==4 % input year only
        yr1=str2double(yr_out);
        yr2=str2double(yr_out);
    elseif length(yr_out)==9 && yr_out(5)=='-' % year range
        yr1=str2double(yr_out(1:4));
        yr2=str2double(yr_out(6:9));
    end
end

% generate time array
t1 = datetime(yr1,1,1,0,0,0);
t2 = datetime(yr2,12,31,23,0,0);
time_array = [t1:hours(1):t2]';

%%% select location and specify time period
if strcmp(loc_in,'EWS') % general DMPs for Eureka
    
    location.latitude=79.989;
    location.longitude=-85.937;
    location.altitude=10;

    % take winter/spring only
    time_array(time_array.Month<10 & time_array.Month>5)=[];
    
elseif strcmp(loc_in,'RL') %DMPs for SOLID measurement times only

    location.latitude=80.053;
    location.longitude=-86.416;
    location.altitude=610;
    
% 24 Jan 2017 - 10 March 2017.
% 18 Oct 2017 - 9 Nov 2017.
% 7 Feb 2018 - 10 March 2018.
% 21 Oct 2018 - 7 Nov 2018.
% % % 3 Feb 2019 - 24 Feb 2019.
% % % 10 Nov 2019 - 3 Dec 2019.
% % %  2 Feb 2020 - 10 March 2020.   
    % include 1 day cushion, just in case
    time_array=time_array((time_array>=datetime(2017,1,23) & time_array<datetime(2017,3,12)) | ...
                          (time_array>=datetime(2017,10,17) & time_array<datetime(2017,11,11)) | ...
                          (time_array>=datetime(2018,2,6) & time_array<datetime(2018,3,12)) | ...
                          (time_array>=datetime(2018,10,20) & time_array<datetime(2018,11,9)) | ...
                          (time_array>=datetime(2019,2,2) & time_array<datetime(2019,2,26)) | ...
                          (time_array>=datetime(2019,10,9) & time_array<datetime(2019,12,5)) | ...
                          (time_array>=datetime(2020,2,1) & time_array<datetime(2020,3,12))    );
    
    % error if years outside SOLID dates
    if isempty(time_array), error(['No solid measurements for ' num2str(yr1) '-' num2str(yr2)]), end
    
end


%%% get time structure
time=struct;

time.year=year(time_array);
time.month=month(time_array);
time.day=day(time_array);
time.hour=hour(time_array);
time.min=minute(time_array);
time.sec=second(time_array);

time.UTC=zeros(size(time.year));


% initialize output structure
LOSinfo=struct;
LOSinfo.z=z;
 
% assign fixed coordinates
% required shape is len(z) x len(time)
LOSinfo.Lat=repmat(location.latitude,length(z),length(time.year));
LOSinfo.Lon=repmat(location.longitude,length(z),length(time.year));



end

