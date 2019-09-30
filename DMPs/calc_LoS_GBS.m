function LOSinfo = calc_LoS_GBS( tag, eu, sza, saa, z_des, sza_array )
%CALC_LOS_GBS Calculates line of sight (LoS) information for GBS instruments
%   Uses fixed list of SZA values
%   INPUT:
%       mode
%           1: LoS info specific to SZA (calculate approx. 
%              LoS info for each SZA-SAA pair)
%       tag ('NO2_UV', 'NO2_VIS', or 'O3_VIS'): select wavelength-dependent
%          scattering heights for each species
%       eu: structure containing latitude, longitude and altitude of
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

% extrapolate scattering heights if list of deisred SZA is passed as input
if nargin==6
    disp('Using extended SZA range')
    
    % create column vector
    if size(sza_array,1)==1, sza_array=sza_array'; end
    
    % assume constant scattering height for SZA < lowest available
    % scattering heigh SZA
    if sza_des(end)==sza_array(end)
        z_scat=[ones(length(sza_array)-length(z_scat),1)*z_scat(1) ; z_scat];
        sza_des=sza_array;
    else
        error('Code up part that handles sza>90')
    end
end

% radius of the Earth in km
R_e = 6378.1;


%% Calculate distance of each LoS point from instrument (in degrees of arc)
% range: rows for each SZA, columns for each altitude 

for i = 1:length(sza_des)
    for j = 1:length(z_des)
        if z_des(j) < z_scat(i)
            range(i,j) = 0;
        else
            range(i,j) = sza_des(i) - asind( (R_e+z_scat(i))./(R_e+z_des(j)) .* sind(sza_des(i)));
        end
    end
end


%% Loop over each measurement

n=0;
for i=1:length(sza)

    % display progress info
    disp_str=['Calculating effective coordinates for datapoint ' num2str(i)];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    

    % find sza closest to measurement sza
    [~,ind_z_scat]=min(abs(sza_des-sza(i)));

    % get effective coordinates
    for j = 1:length(range(1,:))

        % coordinate of point on a sphere given distance from instrument
        % (in degree of arc) and direction (given by saa) 
        [eff_lat, eff_lon] = reckon(eu.latitude, eu.longitude, range(ind_z_scat,j), saa(i));

        % assign values to output structure
        LOSinfo.Lat(j,i)=eff_lat;
        LOSinfo.Lon(j,i)=eff_lon;

    end

end
fprintf('\n')

end

