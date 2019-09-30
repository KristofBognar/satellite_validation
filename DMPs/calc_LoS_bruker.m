function LOSinfo = calc_LoS_bruker(z,SZA,location,alltime)
% This program will output the latitude and longitude of line-of-sight
% altitudes of FTIR measurements using the SZA
%
% INPUTS
% z: must be a vector with the altitude grid desired
% SZA: a vector with the solar zenith angle of the measurements
% location: is a structure with fields for latitude, longitude and altitude
% of the measurement site/instrument location
% alltime: is a structure with fields for year, month, day, hour, minutes
% (min), seconds (sec), and the 
% UTC offset (UTC) for the timezone
%
% OUPUTS
% LOSinfo.Latitude
% LOSinfo.Longitude
% LOSinfo.z
% 
% Adapted from IDL code written by Felicia Kolonjari (2007)
%
% Dan Weaver, University of Toronto
% dweaver@atmosp.physics.utoronto.ca
% November 2015

%% SET UP
%UTCoffset = -5; % ******* THIS SHOULD BE CHANGED FOR NON-EUREKA SITES! ****
Lat1 = location.latitude; % starting location (e.g. where instrument is)
Lon1 = location.longitude; % starting location (e.g. where instrument is)
% Note: longitude W is defined as negative
%
Lat1r = Lat1*pi/180;   			% latitude of site in radians
Long1r = Lon1 * pi/180;			% longitude of site in radians
OC = 6367.45;					% radius of Earth [km]
% 
SZAr = SZA * (pi/180);          % SZA in radians
Nmeas= length(SZA);             % number of measurements 
Nz = length(z);                 % number of altitudes in the grid
% initialize variables
ml = zeros(Nmeas,Nz);mq = zeros(Nmeas,Nz);
AC = zeros(Nmeas,Nz);
Lat2r = zeros(Nmeas,Nz); Lat2 = zeros(Nmeas,Nz);
Long2 = zeros(Nmeas,Nz); Long2r = zeros(Nmeas,Nz);
%%
for i = 1:Nmeas
    time.year = alltime.year(i);
    time.month = alltime.month(i);
    time.day = alltime.day(i);
    time.hour = alltime.hour(i);
    time.min = alltime.min(i);
    time.sec = alltime.sec(i);
    time.UTC = alltime.UTC(i);
    %
    sun = sun_position_mod(time,location);
    sunSZA = sun.zenith; % for verifying a reasonable calculation
    % since the sun_position code returns the angle from the last solar
    % noon, use the end of each day to calculate that day's SN
    time.hour = 24; time.min = 0; time.sec = 0;
%     time.UTC = UTCoffset; 
    %
    sun = sun_position_mod(time,location);
    SolarNoon = 24 - (double(sun.local_hour))/15; % in UTC
    % 
    msmt_time = datenum(0,0,0,alltime.hour(i),alltime.min(i),alltime.sec(i))*24;
    %%
	HA(i) = (msmt_time - SolarNoon) * 15; % multiplying by 15 provides HA in degrees
	HAr(i) = HA(i)* pi/180;
    %%   
    Z = fliplr(z);
	for k = 1:Nz 
		ml = OC/(OC + Z(k));
		mq = sin(pi - SZAr(i)) * ml;
		AC = sin( SZAr(i) - asin( mq ));
        %
		Lat2r(i,k) = asin(sin(Lat1r)*cos(AC) - cos(Lat1r)*sin(AC)*cos(HAr(i)));
        %
		Long2r(i,k) = Long1r - atan2(sin(HAr(i))*sin(AC)*cos(Lat1r), cos(AC)-sin(Lat1r)*sin(Lat2r(i,k)));
        %
    end 
end 
% back to degrees from radians
Lat2 = Lat2r * 180/pi;
Long2 = Long2r * 180/pi;
%%
LOSinfo.z = Z';
LOSinfo.Lat = Lat2';
LOSinfo.Lon = Long2';
%
end

