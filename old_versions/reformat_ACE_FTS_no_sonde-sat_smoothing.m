function reformat_ACE_FTS( tg )
%reformat_ACE_FTS( data, tg ) Reformat ACE-FTS data similar to GBS data format
%   Detailed explanation goes here

%debug
% tg='O3';

AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';

% integration limits, in km
switch tg
    case {'O3','O3UV'}
        lowlim=14;
        highlim=52;
    case 'NO2'
        lowlim=12;
        highlim=40;
end

comment=['_' num2str(lowlim) '-' num2str(highlim) 'km'];

du=2.687e16;

%% Read in data

% location of ACE data -- .mat files will be saved here
datadir = '/home/kristof/work/satellite_validation/ACE-FTS_data/';

% % % Niall's code to read netCDF data, v3p6 filename is hardcoded
% % read_ace_ncdata_for_mat({tg},datadir);

% load data
load([datadir 'ACE_v3p6_' tg '.mat']);


%% filter by geolocation
% only include measurements within 500 km of PEARL (80.05 N, 86.42 W)

range_km=dist_to_PEARL(tanstruct.lat_tangent,tanstruct.lon_tangent);

% filter by distance
goodind=find(range_km<=500);


%% Reformat time data

% convert from MJD to date
mjd2k=tanstruct.date_mjd(goodind)-date2mjd(2000,1,1,0,0,0);
[year, month, day, hour, minute, second] = mjd2date(tanstruct.date_mjd(goodind));

% day of year and fractional day (still starts with 1 on jan 1)
fd=dayofyear(year, month, day, hour, minute, second);

doy=floor(fd);

% am/pm data is in tanstruct.sr1ss0, opposite of GBS ampm mask
ampm=not(logical(tanstruct.sr1ss0(goodind)));

%% Integrate FTS profiles 

% replace fill values with NaN
tanstruct.vmr(tanstruct.vmr==-999)=NaN;
tanstruct.vmr(tanstruct.vmr==0)=NaN;

% air number density (convert to molec/cm3)
num_dens=((6.022e23*tanstruct.pressure_hPa.*100)./(8.314*tanstruct.temperature_K))*1e-6;

% convert VMR to molec/cm3 and only keep filtered values
tmp=tanstruct.vmr.*num_dens;
concentration=tmp(:,goodind);

% error in concentrations
tmp=tanstruct.vmr_error.*num_dens;
error=tmp(:,goodind);


% Output table
ace_fts=table;
ace_fts.year=year';
ace_fts.day=doy';
ace_fts.ampm=ampm';
ace_fts.fd=fd';
ace_fts.mjd2k=mjd2k';
ace_fts.fractional_time=ace_fts.fd-1;
ace_fts.dist=range_km(goodind)';


switch tg
    case 'O3'
        
        % partial column from FTS data
        % only 0.05% diffrence between midpoint and trapezoidal integration
        part_col=integrate(tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint');
        
        %%% AVK smoothing
        
        % DOAS
        load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
        %%% copy paste from VCD code -- TIME IS OFF for VCD code!!!!
        % time when we have the sonde measurement
        sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
        sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)
        % sat measurement time
        measurement_time = ace_fts.year + ace_fts.fractional_time./daysinyear(ace_fts.year);
        sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);
        
        [lut_prof,lut_avk]=read_DOAS_prof_avk(1,[ace_fts.year,ace_fts.day,sonde_in],...
                                              AVK_LUT_dir );
        
        % get avk smoothed partial column
        part_col_smooth=integrate_smooth( tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint',...
                           [0.5:59.5]*1e5, lut_prof(:,4:end), lut_avk(:,4:end), 'interp');      

        % bruker
        [bk_alt,bk_layer_h,bk_ap,bk_avk]=read_bruker_prof_avk(1,ace_fts.mjd2k);

        % get avk smoothed partial column
        % change lowlim to match bruker grid bound nearest to 14 km
        part_col_smooth_bk=integrate_smooth( tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint',...
                           bk_alt*1e5, bk_ap, bk_avk, 'interp',bk_layer_h*1e5);      
    
        % surface to 14km partial column from sonde data
        tot_col=NaN(size(part_col));
        tot_col_smooth=NaN(size(part_col_smooth));
        tot_col_smooth_bk=NaN(size(part_col_smooth_bk));
        col_below_rl=NaN(size(part_col));
        
        for i=1:length(part_col)

            [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                interp_ozonesonde( year(i), fd(i)-1 );

            sonde_conc=sonde_vmr.*sonde_num_dens;
            
            tmp=integrate(sonde_alt_grid*1e2,sonde_conc,1000,61000,'midpoint');
            col_below_rl(i)=tmp;

            tmp=integrate(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint');
            tot_col(i)=part_col(i) + tmp;

            tmp=integrate_smooth(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint',...
                    [0.5:59.5]*1e5, lut_prof(i,4:end), lut_avk(i,4:end), 'layer_mean');
            tot_col_smooth(i)=part_col_smooth(i) + tmp;

            tmp=integrate_smooth(sonde_alt_grid*1e2,sonde_conc,61000,lowlim*1e5,'midpoint',...
                    bk_alt*1e5, bk_ap(i,:), bk_avk(i,:), 'layer_mean',bk_layer_h*1e5);
            tot_col_smooth_bk(i)=part_col_smooth_bk(i) + tmp;
            
            
            
        end
        
        % estimate error 
        part_col_err=sqrt(sum(error(lowlim+1:highlim,:).^2))*1e5;
        
        % update utput table
        ace_fts.part_col=part_col'/du;
        ace_fts.tot_col=tot_col'/du;
        ace_fts.part_col_err=part_col_err'/du;
        ace_fts.tot_col_err=part_col_err'/du; % add error for sonde?? likely very small

        ace_fts.part_col_smooth=part_col_smooth'/du;
        ace_fts.tot_col_smooth=tot_col_smooth'/du;
        ace_fts.part_col_smooth_bk=part_col_smooth_bk'/du;
        ace_fts.tot_col_smooth_bk=tot_col_smooth_bk'/du;
        
        ace_fts.col_below_rl=col_below_rl'/du;
        
    case 'NO2'
        
        part_col=integrate(tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint');
        
        %%% AVK smoothing
        [lut_prof,lut_avk]=read_DOAS_prof_avk(2,[ace_fts.year,ace_fts.fractional_time],...
                                              AVK_LUT_dir );
        
        % get avk smoothed partial column
        part_col_smooth=integrate_smooth(tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint',...
                           [0.5:59.5]*1e5, lut_prof(:,3:end), lut_avk(:,3:end), 'interp');
        % same for UV
        [lut_prof,lut_avk]=read_DOAS_prof_avk(3,[ace_fts.year,ace_fts.fractional_time],...
                                              AVK_LUT_dir );
        
        % get avk smoothed partial column
        part_col_smooth_uv=integrate_smooth(tanstruct.altitude_km*1e5,concentration,...
                           lowlim*1e5,highlim*1e5,'midpoint',...
                           [0.5:59.5]*1e5, lut_prof(:,3:end), lut_avk(:,3:end), 'interp');      

        % estimate error
        part_col_err=sqrt(sum(error(lowlim+1:highlim,:).^2))*1e5;
        
        % update utput table
        ace_fts.part_col=part_col';
        ace_fts.tot_col=part_col';
        ace_fts.part_col_err=part_col_err';
        ace_fts.tot_col_err=part_col_err';
        
        ace_fts.tot_col_smooth=part_col_smooth';
        ace_fts.tot_col_smooth_uv=part_col_smooth_uv';
        
end

% add coordinates
ace_fts.lat=tanstruct.lat_tangent(goodind)';
ace_fts.lon=tanstruct.lon_tangent(goodind)';

save([datadir 'ACE_v3p6_' tg '_table' comment '.mat'],'ace_fts')

end

