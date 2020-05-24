function [time,location,sza,z] = prepDMPinput_brewer( bw_num, yr_out )
%PREPDMPINPUT_BREWER Summary of this function goes here
%   Detailed explanation goes here

%% output variables
time=struct;

% possible locations: Ridge lab abd EWS
loc_rl.latitude=80.053;
loc_rl.longitude=-86.416;
loc_rl.altitude=610;

loc_ews.latitude=79.989;
loc_ews.longitude=-85.937;
loc_ews.altitude=10;

% possible altitude grids
% use same grid as GBSs, with lowest datapoint adapted to instrument location
z_rl=fliplr([0.795, [1.5:1:59.5]]);
z_ews=fliplr([0.5:1:59.5]);


%% specify location based on brewer number and year
loc_changed=0;

if bw_num==69
    
    % Brewer69 location changed over the years!! From Xiaoyi: 'This
    % instrument was in EWS before 26-Feb-2009, and moved to PEARL after
    % that. It has been moved by to EWS since 10-Jul-2013.'
    
    loc_changed=1;
    
    % start and end dates of EWS stay, one row for each period
    % ridge lab data is taken as all data not from EWS
    at_ews=[ [datetime(2001,01,01,00,00,00),datetime(2009,02,25,00,00,00)];...
             [datetime(2013,07,10,00,00,00),datetime(2100,01,01,00,00,00)]    ]; 
    
else
    
    %%% To add other brewers with fixed location: set 'location' and 'z' to
    %%% appropriate location in this bloc
    
    %%% To add other brewers with variable location: follow brewer69 example
    
end

%% load file with brewer data
load(['/home/kristof/work/brewer/brewer_' num2str(bw_num) '_all.mat'])

time_array=brewer_ds.DateTime;

%% break up data
% check if only subset of years is needed
% for fix location, only need to select years
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

% if location has changed, need to restrict index to one location
if loc_changed % variable location: restrict results to one location, give warning

    % data points from EWS
    ind_ews=false(size(time_array));
    for i=1:size(at_ews,1)
        ind_ews=(ind_ews | (brewer_ds.DateTime>at_ews(i,1) & brewer_ds.DateTime<at_ews(i,2)));
    end
    
    % when not at EWS, instrument is assumed to be at the ridge lab
    ind_rl=~ind_ews;
    
    % if data from one location only: proceed as normal
    if sum(ind & ind_ews)==0 % no EWS data
        location=loc_rl;
        z=z_rl;
    elseif sum(ind & ind_rl)==0 % no ridge lab data
        location=loc_ews;
        z=z_ews;

    % if two locations: find which location has more datapoints
    % (considering only the years requested), and process only those data
    elseif sum(ind & ind_ews) >= sum(ind & ind_rl)

        ind=(ind & ind_ews);
        location=loc_ews;
        z=z_ews;

        warning('Multiple instrument locations: only processing data from EWS')
        disp('Periods processed (rows indicate start and end):')
        disp(at_ews)

    elseif sum(ind & ind_ews) < sum(ind & ind_rl)

        ind=(ind & ind_rl);
        location=loc_rl;
        z=z_rl;

        warning('Multiple instrument locations: only processing data from the Ridge Lab')
        disp('Periods shown here were EXCLUDED (rows indicate start and end):')
        disp(at_ews)
    end
end

%% output
% time output
time.year=year(time_array(ind));
time.month=month(time_array(ind));
time.day=day(time_array(ind));
time.hour=hour(time_array(ind));
time.min=minute(time_array(ind));
time.sec=second(time_array(ind));

time.UTC=zeros(size(time.year));

% solar position output
sza=brewer_ds.ZA(ind);

end

