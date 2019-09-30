function [ time, location, sza_out, saa_out, z ] = sza_saa_for_LoS( startyear, endyear, sza_array )
%SZA_SAA_FOR_LOS Calculates time of SZA in sza_array for each day in
%specified time period

%% setup

% default sza where scattering heights are available: include a range of
% SZA, and closest scattering height SZA are found by calc_LoS_GBS.m
if nargin==2
    sza_array=[72:2:90];
end

% timestep for SZA-SAA pairs in minutes
timestep=10;

% variables
time=struct;
sza=[];
saa=[];

% location info
location.latitude=80.053;
location.longitude=-86.416;
location.altitude=610;

z=[0.795, [1.5:1:59.5]];
if z(1) < z(length(z))
    z = fliplr(z);
end

%% calculate solar positin for all the years

disp('Calculating solar position for specified period');

% UTC time interval for SZA calculation (start with local midnight)
t1=datetime(startyear,02,01,06,00,00);
t2=datetime(endyear,11,01,05,60-timestep,00);

time_array=t1:minutes(timestep):t2;
time_array=time_array';


% get SZA and SAA for entire time interval
[saa,el]=SolarAzEl(time_array,...
                  repmat(location.latitude, length(time_array),1),...
                  repmat(location.longitude, length(time_array),1),...
                  repmat(location.altitude/1000, length(time_array),1));

% convert elevation to SZA
sza=90-el;


%% loop over days to pick out relevant SZA

disp('Calculating time of specified SZA for each day');

loop_step=24*60/timestep; % loop to move one day at a time
ind_keep=[]; % indices of sza array to keep
for i=1:loop_step:length(time_array)-1
    
    % indices for given day
    ind=i:i+loop_step-1;
    
    % range in sza
    % leave a 1 deg margin since SZA step is 2deg, so up to 1 deg deviation is fine
    sza_min=min(sza(ind))-1;
    sza_max=max(sza(ind))+1;
    
    % skip early spring/late fall
    if sza_min>max(sza_array), continue, end
    
    % only consider spring prior to 2005 (when SAOZ was installed)
    % UT-GBS only started year-round measurements in 2008
    if year(time_array(ind(1)))<2005 && month(time_array(ind(1)))>4, continue, end
    
    % loop over scattering height SZA
    for j=sza_array

        % check if sza exists on given day
        if j>=sza_min && j<=sza_max
            
            % find times for sza (no interpolation, assume time resolution
            % is good enough to pick 2 closest times -- given by first 2
            % indices of sortind)
            tmp=ones(size(sza(ind)))*j;
            [~,sortind]=sort(abs(sza(ind)-tmp));
            
            % adjust for indexing indexed variables, and 90<sza<91 (only
            % need noon time and saa)
            if sza_min>=89,
                tmp=ind(1)+sortind(1)-1;
            else
                tmp=ind(1)+sortind(1:2)-1;
            end
            
            % indices to keep
            ind_keep=[ind_keep; tmp];
    
        end
    end
end


%% output

% sort by time
[~,timeind]=sort(time_array(ind_keep));
ind_keep=ind_keep(timeind);

% assign final values
sza_out=sza(ind_keep);
saa_out=saa(ind_keep);

time.year=year(time_array(ind_keep));
time.month=month(time_array(ind_keep));
time.day=day(time_array(ind_keep));
time.hour=hour(time_array(ind_keep));
time.min=minute(time_array(ind_keep));
time.sec=second(time_array(ind_keep));

time.UTC=zeros(size(time.year));

end

