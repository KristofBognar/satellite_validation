function [ LOSinfo, time_out ] = calc_LoS_GBS_general( tag, eu, z_des, years )
%CALC_LOS_GBS_GENERAL Summary of this function goes here
%   Detailed explanation goes here


error('Problem with time calculation -- fix it first');

%% Setup

% initialize output structure
LOSinfo=struct;
LOSinfo.z=z_des;

time=struct;
time_out = struct;

% select scattering heights and corresponding SZA (adapted from Cristen's
% thesis and from calc_eff_coord2.m)
switch tag
    case 'NO2_UV'
        scat_vec= [86,16.1;88,19.1;90,23.5];
    case 'NO2_VIS'
        scat_vec= [86,12.4;88,15.1;90,19.5];
    case 'O3_VIS'
        scat_vec= [86,9.8;88,11.7;90,15.6];
    otherwise
        error('First input can only be NO2_UV, NO2_VIS, or O3_VIS')
end

% assign SZA and scattering height
sza_des = scat_vec(:,1);
z_scat = scat_vec(:,2);

% radius of the Earth in km
R_e = 6378.1;

day_beg=1;
day_end=365;

% % % % %% time interval for SZA calculation
% % % % t1=datetime(years(1),01,01,00,00,00);
% % % % t2=datetime(years(2),01,01,00,00,00);
% % % % 
% % % % time_array=t1:minutes(10):t2;
% % % % time_array=time_array';
% % % % 
% % % % 
% % % % %% get SZA and SAA for entire time interval
% % % % [az,el]=SolarAzEl(time_array,...
% % % %                   repmat(eu.latitude, length(time_array),1),...
% % % %                   repmat(eu.longitude, length(time_array),1),...
% % % %                   repmat(eu.altitude, length(time_array),1));
% % % % 
% % % % 
% % % % sza_list=[az,90-el];


%% loop over years
count=0;
for yyyy=years(1):years(2)
    
    time.year=yyyy;
    time.UTC = 0;
    time.sec = 0;


    % build SZA versus SAA table for this time-period
%     disp('Building SZA/azimuth lookup table')
    sza_list = [];
    for jjj = day_beg:day_end,
        for hour = 0:23,
            for min = 0:10:50,
                [day, month] = Julian2Date(yyyy, jjj);
                time.day = day;
                time.month = month;
                time.hour = hour;
                time.min = min;
                param = sun_position(time,eu);
                fd = jjj + time.hour / 24 + time.min / 60;
                sza_list = [sza_list; jjj fd param.azimuth param.zenith];
            end
        end
    end

    % Now sort out the local day issue
    for jjj = day_beg:day_end
        i_day = find(sza_list(:,1) == jjj);
        sza_day = sza_list(i_day,:);
        [~, i_midnight] = max(sza_day(:, 3)); 
        sza_day(1:i_midnight,1) = (sza_day(1,1) - 1)*ones(i_midnight,1);
        sza_list(i_day,:) = sza_day;
    end

    % now with our fabulous list make calculations
    sza_list_i = [];
    for jjj = day_beg:day_end
        for ampm = 0:1,
             % get morning or evening indices for the given day
            if ampm == 0
                ind = find(sza_list(:,1) == jjj & sza_list(:,3) < 180);
            elseif ampm == 1
                ind = find(sza_list(:,1) == jjj & sza_list(:,3) > 180);
            end
            for sza_i = sza_des'
                % grab the relevant data and interpolate to the sza grid
                sza_list_i = [sza_list_i;
                    [jjj ampm...
                    interp1(sza_list(ind,4), [sza_list(ind,2:3)], sza_i)...
                    sza_i]];
            end
        end
    end

    % so now we have a lovely list of values interpolated to the desired sza!
    % let's reorganize
    ind = find(~isnan(sza_list_i(:,4)));
    sza_list = sza_list_i(ind,:);

    % calculate the range
    for i = 1:length(sza_des)
        for j = 1:length(z_des)
            if z_des(j) < z_scat(i)
                range(i,j) = 0;
            else
                range(i,j) = sza_des(i) - asind( (R_e+z_scat(i))./(R_e+z_des(j)) .* sind(sza_des(i)));
            end
        end
    end

    % and calculate the longitude and latitude
    for i = 1:length(sza_list(:,1))

        % get dates and times
        [dd, mm] = Julian2Date(yyyy, sza_list(i,1));
        time_ft = (sza_list(i,3) - floor(sza_list(i,3)))*24*60*60;
        time_ft = floor(time_ft);

        for j = 1:length(range(1,:))
            k = find(sza_des == sza_list(i,5));
            if isnan(range(k,j)), continue; end
            [eff_lat, eff_lon] = reckon(eu.latitude, eu.longitude, range(k,j), sza_list(i,4));

            % assign values to output structure
            LOSinfo.Lat(j,count+i)=eff_lat;
            LOSinfo.Lon(j,count+i)=eff_lon;
            
            time_out.year(count+i)=yyyy;
            time_out.month(count+i)=mm;
            time_out.day(count+i)=dd;
            time_out.hour(count+i)=floor(time_ft/3600);
            time_out.min(count+i)=floor(((time_ft/3600)-time_out.hour(count+i))*60);
            time_out.sec(count+i)=0;
            
        end
    end

%     range_km = range * (pi/180) * R_e;
    
    count=count+i;

end

end

