function LOSinfo = calc_LoS_GBS_interp( tag, location, sza, saa, z_des )
%CALC_LOS_GBS Calculates line of sight (LoS) information for GBS instruments
%   Uses solar position of each measurement, scattering heights are interpolated
%   INPUT:
%       tag ('NO2_UV', 'NO2_VIS', or 'O3_VIS'): select wavelength-dependent
%          scattering heights for each species
%       location: structure containing latitude, longitude and altitude of
%          instrument
%       sza: solar zenith angle of the measurements
%       saa: solar azimuth angle of the measurements
%       z_des: altitude grid for LoS caculation
%
%   OUTPUT:
%       LOSinfo.Lat: latitude of point with altitude z along LoS (length(z) x length(sza))
%       LOSinfo.Lon: longitude of point with altitude z along LoS (length(z) x length(sza))
%       LOSinfo.z: returns z_des

%% Setup

% radius of the Earth in km
R_e = 6378.1;

% initialize output structure
LOSinfo=struct;
LOSinfo.z=z_des;

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


% assume constant scattering height for SZA < 86 and SZA > 90
sza_tmp=sza;
sza_tmp(sza_tmp>90)=90;
sza_tmp(sza_tmp<86)=86;

% interpolate scattering heights for each measurement
z_scat_interp=interp1(sza_des,z_scat,sza_tmp);
    

%% Loop over measurements
n=0;
for i=1:length(sza)

    % display progress info
    disp_str=['Calculating effective coordinates for datapoint ' num2str(i)];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    

    % loop over altitudes
    for j = 1:length(z_des)
        %% Calculate distance of each LoS point from instrument (in degrees of arc)
        % range: rows for each SZA, columns for each altitude 
        
        if z_des(j) < z_scat_interp(i)
            range(j,i) = 0;
        else
            range(j,i) = sza(i) - asind( (R_e+z_scat_interp(i))./(R_e+z_des(j)) .* sind(sza(i)));
        end
        
        %% get effective coordinates
        
        % coordinate of point on a sphere given distance from instrument
        % (in degree of arc) and direction (given by saa) 
        [eff_lat, eff_lon] = reckon(location.latitude, location.longitude, range(j,i), saa(i));

        % assign values to output structure
        LOSinfo.Lat(j,i)=eff_lat;
        LOSinfo.Lon(j,i)=eff_lon;

    end

end
fprintf('\n')

end

