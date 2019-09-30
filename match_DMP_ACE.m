function [ spv, temperature, theta, lat, lon ] = match_DMP_ACE( mjd2k, alt )
%MATCH_DMP_DOAS Returns DMPs for given DOAS measurement
%
% INPUT
%       mjd2k_in: array of measurement times ( Jan, 1, 2000, 00:00 = 0)
%       alt_in: altitude levels where DMP outputs are required
%       !! dmp altitude grid is [0.5,1.5,2.5,...,74.5]
%
% OUTPUT
%       spv, temperature, theta: DMPs for the given measurement,
%               dimensions of length(ft_meas) X length(alt)
%       lat, lon: latitude and longitude of the point along the line of
%               sight at each altitude, dimensions of length(ft_meas) X length(alt)
%
% Kristof Bognar, May 2018

if size(alt,1)~=1, alt=alt'; end

%% setup output variables

spv=NaN(length(mjd2k),length(alt));
temperature=NaN(length(mjd2k),length(alt));
lat=NaN(length(mjd2k),length(alt));
lon=NaN(length(mjd2k),length(alt));
theta=NaN(length(mjd2k),length(alt));

%% pair to the measurements

load('/home/kristof/work/DMP/DMP_ACE_failed_GLC_filled/ACE_DMP_struct.mat');
    
dmp_date=dmpstruct.mjd2k;

% structure using Niall's code - cannot read netCDF DMPs
% load('/home/kristof/work/DMP/ACE_v3p5_6_DMP2p0_all.mat');
% dmp_date=dmpstruct.date_mjd-51544;

for i=1:length(mjd2k)

    % find index of closest time
    [tmp,ind_dmp]=sort(abs(dmp_date-mjd2k(i)));
    
    % check if it's an actual match (actual diff sould be <1e-5, i.e. <1 second)
    if tmp(1)<1e-4 
        
        % index of match
        ind_dmp=ind_dmp(1);

        spv(i,:)=interp1(dmpstruct.altitude_km(:,ind_dmp),dmpstruct.spv(:,ind_dmp),alt);
        temperature(i,:)=interp1(dmpstruct.altitude_km(:,ind_dmp),dmpstruct.T(:,ind_dmp),alt);
        lat(i,:)=interp1(dmpstruct.altitude_km(:,ind_dmp),dmpstruct.lat(:,ind_dmp),alt);
        lon(i,:)=interp1(dmpstruct.altitude_km(:,ind_dmp),dmpstruct.lon(:,ind_dmp),alt);
        theta(i,:)=interp1(dmpstruct.altitude_km(:,ind_dmp),dmpstruct.Theta(:,ind_dmp),alt);

    end
    
end

