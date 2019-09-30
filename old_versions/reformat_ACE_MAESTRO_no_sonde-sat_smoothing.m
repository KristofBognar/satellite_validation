function reformat_ACE_MAESTRO( tg )
%REFORMAT_ACE_MAESTRO(tg) reformat ACE-MAESTRO data as a matlab table; intgrate columns 
%
% Matches MAESTRO data to FTS measurements and uses FTS P, T data to
% calculate number densities
%
% Data filters:
%   Distance from PEARL (<=500 km)
%   No orbit match in ACE-FTS data (noted if measurement was within 500 km of PEARL)
%
% Unfiltered (to have some idea about the % of bad measurements near PEARL):
%   Bad values (returns NaN column)
%   Layers where retrieval failed (cout between altitude limits returned)


% integration limits, in km
switch tg
    case {'O3','O3UV'}
        lowlim=14;
        highlim=52;
    case 'NO2'
        lowlim=12;
%         lowlim=17;
        highlim=40;
end

% add extra infor to filename
% comment='17-40km';
comment=['_' num2str(lowlim) '-' num2str(highlim) 'km'];

du=2.687e16;

data_dir='/home/kristof/work/satellite_validation/ACE-MAESTRO_data/';
ace_dir = '/home/kristof/work/satellite_validation/ACE-FTS_data/';
cur_dir=pwd;

% load ACE FTS data to match orbits
% orbit, P, T data is the same for O3 and NO2 in FTS data
load([ace_dir 'ACE_v3p6_O3.mat'])

% read geolocation files -- only used to double check orbit match
load([data_dir 'sunrise_table.mat'])
load([data_dir 'sunset_table.mat'])

% change to MAESRTO dir
cd([data_dir tg]);

% Make list of filenames
tmp=dir('*.dat');
flist={tmp.name};

% sr/ss, orbit arrays
sr1ss0=NaN(1,length(flist));
occultation=NaN(1,length(flist));

% define output variables
no_fts_orbit=[];
bad_files=0;

ace_maestro=table;
ace_maestro.year=0;
ace_maestro.day=0;
ace_maestro.ampm=0;
ace_maestro.fd=0;
ace_maestro.mjd2k=0;
ace_maestro.fractional_time=0;
ace_maestro.dist=0;
ace_maestro.part_col=0;
ace_maestro.tot_col=0;
ace_maestro.part_col_err=0;
ace_maestro.tot_col_err=0;
ace_maestro.part_col_smooth=0;
ace_maestro.tot_col_smooth=0;
ace_maestro.part_col_smooth_bk=0;
ace_maestro.tot_col_smooth_bk=0;
ace_maestro.col_below_rl=0;
ace_maestro.lat=0;
ace_maestro.lon=0;
ace_maestro.occultation=0;
ace_maestro.failed_retr_frac=0;

% stuff for AVK smoothing
AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';
load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
%%% copy paste from VCD code -- TIME IS OFF for VCD code!!!!
% time when we have the sonde measurement
sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)

% Loop over files
n=0;
for i=1:length(flist)
    %% display progress info
    disp_str=['Reading file ', num2str(i), '/', num2str(length(flist))];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
    
    %% get orbit and sunrise/sunset info
    tmp=strsplit(char(flist{i}),'_');
    tmp=tmp{1};
    
    % sunrise/sunset mask
    if strcmp(tmp(1:2),'sr')
        sr1ss0(i)=1;
    elseif strcmp(tmp(1:2),'ss')
        sr1ss0(i)=0;
    end
    
    % occultation number
    occultation(i)=str2double(tmp(3:end));

    % find matching orbit in FTS data
    orbit_fts=find(tanstruct.occultation==occultation(i) &...
                   tanstruct.sr1ss0 == sr1ss0(i));
    
    % check if there's multiple matches
    if length(orbit_fts)>1 || isempty(orbit_fts)
        
        % try to select orbit based on maestro look-up table coords
        if sr1ss0(i)==1
            lon_lut=sunrise_table.longitude(sunrise_table.orbit_number==occultation(i));
            lat_lut=sunrise_table.latitude(sunrise_table.orbit_number==occultation(i));
        elseif sr1ss0(i)==0
            lon_lut=sunset_table.longitude(sunset_table.orbit_number==occultation(i));
            lat_lut=sunset_table.latitude(sunset_table.orbit_number==occultation(i));
        end
        
        orbit_fts=find(tanstruct.occultation==occultation(i) &...
                       tanstruct.sr1ss0 == sr1ss0(i) &...
                       tanstruct.lon_tangent-lon_lut < 0.1);
        % give up
        if length(orbit_fts)>1, continue, end
        
        % save orbit if close to pearl but no FTS match was foud
        if isempty(orbit_fts)
            range_tmp=dist_to_PEARL(lat_lut, lon_lut);
            if range_tmp>500
                continue
            else
                no_fts_orbit=[no_fts_orbit;[occultation(i),sr1ss0(i),range_tmp]];
                continue
            end
        end
    end
    
    %% filter by geolocation
    % only include measurements within 500 km of PEARL (80.05 N, 86.42 W)

    range_km=dist_to_PEARL( tanstruct.lat_tangent(orbit_fts),...
                            tanstruct.lon_tangent(orbit_fts)    );

    % filter by distance
    if range_km>500, continue, end

    %% calculate number density profiles
    try
        data=read_maestro(flist{i});
    catch
        bad_files=bad_files+1;
        continue
    end
        
    % interpolate VMR data to FTS grid
    tmp_prof=interp1(data(:,2),data(:,3),tanstruct.altitude_km);
    
    % air number density (convert to molec/cm3)
    num_dens=((6.022e23*tanstruct.pressure_hPa(:,orbit_fts).*100)./...
              (8.314*tanstruct.temperature_K(:,orbit_fts)))*1e-6;
    
    % calculate number density (molec/cm3)
    concentration=tmp_prof.*num_dens;

    % error in concentrations
    tmp_err=interp1(data(:,2),data(:,4),tanstruct.altitude_km);    
    error=tmp_err.*num_dens;
    
    %% Reformat time data

    % convert from MJD to date
    mjd2k=tanstruct.date_mjd(orbit_fts)-date2mjd(2000,1,1,0,0,0);
    [year, month, day, hour, minute, second] = mjd2date(tanstruct.date_mjd(orbit_fts));

    % day of year and fractional day (still starts with 1 on jan 1)
    fd=dayofyear(year, month, day, hour, minute, second);

    doy=floor(fd);

    % am/pm data is in tanstruct.sr1ss0, opposite of GBS ampm mask
    ampm=not(logical(sr1ss0(i)));

    %% save some data
    tmp_table=table;
    tmp_table.mjd2k=mjd2k;
    tmp_table.year=year;
    tmp_table.day=doy;
    tmp_table.ampm=ampm;
    tmp_table.fd=fd;
    tmp_table.fractional_time=fd-1;
    tmp_table.dist=range_km;

    %% integrate profiles
    
    % check if all data in altitude range is retrieved, skip if not
%     if any(data(data(:,2)<=highlim & data(:,2)>=lowlim,5)==0), continue, end
    % get number of layers where data is a-priori
    no_retr=sum(not(logical(data(data(:,2)<=highlim & data(:,2)>=lowlim,5))));
    tot_layers=sum(data(:,2)<=highlim & data(:,2)>=lowlim);

    % partial column from FTS data
    % only 0.05% diffrence between midpoint and trapezoidal integration
    part_col=integrate(tanstruct.altitude_km*1e5,concentration,...
                       lowlim*1e5,highlim*1e5,'midpoint');

    % estimate error 
    part_col_err=sqrt(sum(error(lowlim+1:highlim,:).^2))*1e5;

    % error sign check
    if any(error(lowlim+1:highlim,:)<0),
        bad_error=true;
    else
        bad_error=false;
    end


    switch tg
        case {'O3', 'O3UV'}

            %%% AVK smoothing
            % DOAS
            % sat measurement time
            measurement_time = tmp_table.year + tmp_table.fractional_time./daysinyear(tmp_table.year);
            sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

            [lut_prof,lut_avk]=read_DOAS_prof_avk(1,[tmp_table.year,tmp_table.day,sonde_in],...
                                                  AVK_LUT_dir );

            % get avk smoothed partial column
            part_col_smooth=integrate_smooth( tanstruct.altitude_km*1e5,concentration,...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               [0.5:59.5]*1e5, lut_prof(:,4:end), lut_avk(:,4:end), 'interp');      
            
            % bruker
            [bk_alt,bk_layer_h,bk_ap,bk_avk]=read_bruker_prof_avk(1,tmp_table.mjd2k);

            % get avk smoothed partial column
            % change lowlim to match bruker grid bound nearest to 14 km
            part_col_smooth_bk=integrate_smooth( tanstruct.altitude_km*1e5,concentration,...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               bk_alt*1e5, bk_ap, bk_avk, 'interp',bk_layer_h*1e5);      
                           
            % surface to 14km partial column from sonde data
            if part_col>0 
                [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                    interp_ozonesonde( year, fd-1 );

                sonde_conc=sonde_vmr.*sonde_num_dens;

                tmp=integrate(sonde_alt_grid*1e2,sonde_conc,1000,61000,'midpoint');
                col_below_rl=tmp;

                tmp=integrate(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint');
                tot_col=part_col + tmp;

            else % bad value placeholder (-1.#JE+000) read as -1 by code
                part_col=NaN;
                tot_col=NaN;
                part_col_err=NaN;
                col_below_rl=NaN;
            end
            
            % surface to 14km partial column from sonde data
            % do it separately , in case smoothing changes partial column
            if part_col_smooth>0 
                [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                    interp_ozonesonde( year, fd-1 );

                sonde_conc=sonde_vmr.*sonde_num_dens;

                tmp=integrate_smooth(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint',...
                    [0.5:59.5]*1e5, lut_prof(4:end), lut_avk(4:end), 'layer_mean');
                
                tot_col_smooth=part_col_smooth + tmp;

            else % bad value placeholder (-1.#JE+000) read as -1 by code
                part_col_smooth=NaN;
                tot_col_smooth=NaN;
            end
            
            % surface to 14km partial column from sonde data
            % do it separately , in case smoothing changes partial column
            if part_col_smooth_bk>0 
                [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                    interp_ozonesonde( year, fd-1 );

                sonde_conc=sonde_vmr.*sonde_num_dens;

                tmp=integrate_smooth(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint',...
                    bk_alt*1e5, bk_ap, bk_avk, 'layer_mean',bk_layer_h*1e5);
                
                tot_col_smooth_bk=part_col_smooth_bk + tmp;

            else % bad value placeholder (-1.#JE+000) read as -1 by code
                part_col_smooth_bk=NaN;
                tot_col_smooth_bk=NaN;
            end
            
            if bad_error, part_col_err=NaN; end

            % update utput table
            tmp_table.part_col=part_col/du;
            tmp_table.tot_col=tot_col/du;
            tmp_table.part_col_err=part_col_err/du;
            % add error for sonde?? likely very small
            tmp_table.tot_col_err=part_col_err/du;
            
            tmp_table.part_col_smooth=part_col_smooth/du;
            tmp_table.tot_col_smooth=tot_col_smooth/du;
            tmp_table.part_col_smooth_bk=part_col_smooth_bk/du;
            tmp_table.tot_col_smooth_bk=tot_col_smooth_bk/du;
            
            tmp_table.col_below_rl=col_below_rl/du;
            
        case 'NO2'

            if part_col<0 
                part_col=NaN;
                part_col_err=NaN;
            end
            
            if bad_error, part_col_err=NaN; end
            
            % update utput table
            tmp_table.part_col=part_col;
            tmp_table.tot_col=part_col;
            tmp_table.part_col_err=part_col_err;
            tmp_table.tot_col_err=part_col_err;

    end
    
    % add failed retrieval fraction
    tmp_table.failed_retr_frac=no_retr/tot_layers;
    
    % add coordinates
    tmp_table.lat=tanstruct.lat_tangent(orbit_fts);
    tmp_table.lon=tanstruct.lon_tangent(orbit_fts);
    tmp_table.occultation=occultation(i);

    % add to output table
    ace_maestro=[ace_maestro;tmp_table];
    
end

%% save data
% remove line left over from table setup
ace_maestro(1,:)=[];

% sort data by date
ace_maestro=sortrows(ace_maestro,[1 2 3]);

% save
save([data_dir 'ACE_MAESTRO_' tg '_table' comment '.mat'],'ace_maestro')

fprintf('\n')
disp('Done')

if ~isempty(no_fts_orbit)
    disp(['Could not find the following orbits'])
    no_fts_orbit
end

cd(cur_dir);

bad_files

end

