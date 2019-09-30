function [ spv, temperature, theta, lat, lon ] = match_DMP_PARIS( ft_meas, year, alt )
%MATCH_DMP_PARIS Returns DMPs for given PARIS measurement
%
% INPUT
%       ft_meas: array of measurement times (fractional time, Jan, 1, 00:00 = 0)
%       year: year of measurments as a number, or as an array with the same
%               length as ft_meas
%       alt: altitude levels where DMP outputs are required
%       !! dmp altitude grid is [0.8,1.5,2.5,...,59.5]
%
% OUTPUT
%       spv, temperature, theta: DMPs for the given measurement,
%               dimensions of length(ft_meas) X length(alt)
%       lat, lon: latitude and longitude of the point along the line of
%               sight at each altitude, dimensions of length(ft_meas) X length(alt)
%
%       Note: output will match input line by line, even if input times are not sorted
%
% Kristof Bognar, May 2018

% if ~issorted(year), error('measurements must be sorted'); end

if size(alt,1)~=1, alt=alt'; end

%% check how year is specified
if length(year)==1 % sinle year
    
    % set loop variables
    years=year;
    ind_ft=1:length(ft_meas);
    
elseif length(year)==length(ft_meas) % year specified for each fractional time
    
    % set loop variables
    [years,~,ind_ft_unique]=unique(year);
    
    % make sure shape is correct
    if size(years,1)~=1, years=years'; end
    if size(ind_ft_unique,1)~=1, ind_ft_unique=ind_ft_unique'; end
    
else
    error('Measurement times and years don''t match')
end

%% setup output variables

spv=NaN(length(ft_meas),length(alt));
temperature=NaN(length(ft_meas),length(alt));
lat=NaN(length(ft_meas),length(alt));
lon=NaN(length(ft_meas),length(alt));
theta=NaN(length(ft_meas),length(alt));

%% loop over the measurements
n=0;
for yy=years
    
    % display progress info
    disp_str=['DMP matching for ' num2str(yy) ' PARIS data'];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
    
    % select list of indices for given year if there are multiple years
    if exist('ind_ft_unique', 'var')
        ind_ft=find(ind_ft_unique==find(years==yy));
    end
    
    % load DMP file for given year and species
    % saved using /home/kristof/work/matlab_code/DMPs/reformat_DMPs.m
    try
        load(['/home/kristof/work/DMP/DMP_PARIS/PARIS_DMP_table_' num2str(yy) '.mat']);
    catch
        continue
    end
    
    for ii=ind_ft
        
        [tmp,ind_dmp]=sort(abs(fractional_time-ft_meas(ii)));
        
        % check if it's an actual match (actual diff sould be <1e-5, i.e. <1 second)        
        if tmp(1)<1e-4
            
            % index of match
            ind_dmp=ind_dmp(1);

            spv(ii,:)=interp1(dmp_all{ind_dmp}.alt,dmp_all{ind_dmp}.spv,alt);
            temperature(ii,:)=interp1(dmp_all{ind_dmp}.alt,dmp_all{ind_dmp}.temperature,alt);
            lat(ii,:)=interp1(dmp_all{ind_dmp}.alt,dmp_all{ind_dmp}.lat,alt);
            lon(ii,:)=interp1(dmp_all{ind_dmp}.alt,dmp_all{ind_dmp}.lon,alt);
            theta(ii,:)=interp1(dmp_all{ind_dmp}.alt,dmp_all{ind_dmp}.theta,alt);
            
        end
    end
end

fprintf('\n')

end

