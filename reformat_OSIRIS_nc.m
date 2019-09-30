function reformat_OSIRIS_nc( tg )
%REFORMAT_OSIRIS(tg): reformat OSIRIS HDF data as a matlab table


switch tg
    case {'O3','O3UV'}
        lowlim=14;
        highlim=52;
    case 'NO2'
        lowlim=12;
        highlim=32; 
        % mean top alt is 36.58 km; with 40 km limit, 85% of measurements
        % are missing the top; the number is 22.5 and 12 for 36 and 35 km, respectively
        % missing column: 40 km - 1.8%; 36km - 7.6%; 35km - 10% -- more
        % details in notes
end

% max acceptable distance from PEARL in km
comparison_range=500;

disp(['Using ' num2str(comparison_range) ' km radius around PEARL']);

comment=['_' num2str(lowlim) '-' num2str(highlim) 'km'];

% stuff for AVK smoothing
AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';
load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
%%% copy paste from VCD code -- TIME IS OFF for VCD code!!!!
% time when we have the sonde measurement
sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)


%% make list of data files

cur_dir=pwd;
data_dir='/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/';

cd([data_dir tg]);

tmp=dir('*.nc');
flist={tmp.name};
% flist(1:2)=[];

%% setup output table
osiris=table;   


%% loop through files and read in/reformat data

year_prev=2000; % for progress display
du=2.687e16;

top_alt=[];

for i=1:length(flist)
    
    % load variables
    altitude=ncread(flist{i},'altitude'); % in km
    latitude=ncread(flist{i},'latitude');
    longitude=ncread(flist{i},'longitude');
    time=ncread(flist{i},'time');

    % filter by geolocation (only keep measurements within 500km of PEARL)
    range_km=dist_to_PEARL(latitude,longitude);
    goodind=find(range_km<=comparison_range);

    % skip if there are no measurements close to PEARL
    if isempty(goodind), continue, end
    
    % get time info
    % time is days since 1900-01-01 00:00:00
    mjd2k=time(goodind)+ft_to_mjd2k(0,1900);
    dates=mjd2k_to_date(mjd2k);
    
    [fractional_time]=fracdate(dates);
    [year, month, day]=ymd(dates);
        
    % print progress info
    if year_prev~=year(1)
        disp(['Reading ' num2str(year(1)) ' data']);
        year_prev=year(1);
    end
    
    % store results in table
    tmp_table=table;
    tmp_table.mjd2k=mjd2k;
    tmp_table.year=year;
    tmp_table.month=month;
    tmp_table.day=day;
    tmp_table.fractional_time=fractional_time;
    tmp_table.dist=range_km(goodind);
    
    
    %% integrate profiles
    
    switch tg
        case 'O3'
            
            % load ozone data
            tg_num_dens=ncread(flist{i},'ozone_concentration'); % in mole/m3
            tg_error=ncread(flist{i},'ozone_concentration_standard_error'); % in mole/m3

            % convert to molec/cm3
            tg_num_dens=tg_num_dens*6.022140857e+17; 
            tg_error=tg_error*6.022140857e+17;
            
% %             % get P, T
% %             pressure_hPa=ncread(flist{i},'pressure'); 
% %             temperature_K=ncread(flist{i},'temperature'); 
% %             % air number density (convert to molec/cm3)
% %             num_dens=((6.022e23*pressure_hPa.*100)./(8.314*temperature_K))*1e-6;           
% %             % get VMR
% %             vmr=O3_num_dens./num_dens;
            
            % satellite partial column
            part_col=integrate(altitude*1e5, tg_num_dens(:,goodind),...
                               lowlim*1e5, highlim*1e5, 'midpoint');
                           
            % set up array for integration
            concentration=NaN(60,length(goodind));
            concentration(11:end,:)=tg_num_dens(:,goodind);
            
            alt_ext=[[0.5:9.5]';altitude];
            
            %%% AVK smoothing
            % DOAS
            % sat measurement time
            measurement_time = tmp_table.year + tmp_table.fractional_time./daysinyear(tmp_table.year);
            sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

            [lut_prof,lut_avk]=read_DOAS_prof_avk(1,[tmp_table.year,...
                                                  floor(tmp_table.fractional_time)+1,...
                                                  sonde_in], AVK_LUT_dir );

            % bruker
            [bk_alt,bk_layer_h,bk_ap,bk_avk]=read_bruker_prof_avk(1,tmp_table.mjd2k);

            % PARIS
            [pa_alt,pa_layer_h,pa_ap,pa_avk]=read_bruker_prof_avk(1,tmp_table.mjd2k,'PARIS');

            % surface to ridge lab partial column from sonde data
            % subtract from total column (add later for Brewer comparisons)
            col_below_rl=NaN(size(part_col));
            
            % get Bruker avk smoothed partial column 
            part_col_smooth_bk=integrate_smooth( alt_ext*1e5,concentration,...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               bk_alt'*1e5, bk_ap', bk_avk', 'interp',bk_layer_h*1e5);      
            
            % get PARIS avk smoothed partial column 
            part_col_smooth_pa=integrate_smooth( alt_ext*1e5,concentration,...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               pa_alt'*1e5, pa_ap', pa_avk', 'interp',pa_layer_h*1e5);      
            
            
            for j=1:length(part_col)

                [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                    interp_ozonesonde( year(j), fractional_time(j) );

                sonde_conc=sonde_vmr.*sonde_num_dens;

                tmp=integrate(sonde_alt_grid*1e2,sonde_conc,1000,61000,'midpoint');
                col_below_rl(j)=tmp;

                % get effective layer values from sonde, and fill up/replace
                % values in satellite profiles (FTS grid is standard 1km)
                for k=1:lowlim
                    ind_tmp=find(sonde_alt_grid>(k-1)*1000 & sonde_alt_grid<k*1000);
                    concentration(k,j)=nanmean(sonde_conc(ind_tmp));
                end

                % smooth out possible kink with moving average
                boxcar_smoothed=boxcar(alt_ext,concentration(:,j),1);
                % use 2 points on either side of lowlim (14 km)
                concentration(lowlim-1:lowlim+2,j)=boxcar_smoothed(lowlim-1:lowlim+2);
            
            end
            
            % total column from FTS+sonde data
            tot_col=integrate(alt_ext*1e5,concentration,...
                               0,highlim*1e5,'midpoint');

            % get DOAS avk smoothed column (start from RL alt since bruker AVKs start from there)
            tot_col_smooth=integrate_smooth( alt_ext*1e5,concentration,...
                               61000,highlim*1e5,'midpoint',...
                               [0.5:59.5]'*1e5, lut_prof(:,4:end)', lut_avk(:,4:end)', 'interp');      

            % get Bruker avk smoothed column (start from RL alt since bruker AVKs start from there)
            tot_col_smooth_bk=integrate_smooth( alt_ext*1e5,concentration,...
                               61000,highlim*1e5,'midpoint',...
                               bk_alt'*1e5, bk_ap', bk_avk', 'interp',bk_layer_h*1e5);      
            
            % get PARIS avk smoothed column (start from RL alt since bruker AVKs start from there)
            tot_col_smooth_pa=integrate_smooth( alt_ext*1e5,concentration,...
                               62000,highlim*1e5,'midpoint',...
                               pa_alt'*1e5, pa_ap', pa_avk', 'interp',pa_layer_h*1e5);      
            
            % correct for ridge lab altitude
            tot_col=tot_col-col_below_rl;
                % smoothed columns start from RL altitude -- correction done by
                % integration code

            
            % calculate error
            error=tg_error(:,goodind);
              
            % estimate column error (14-52 km limits)
            part_col_err=sqrt(sum(error(lowlim-9:highlim-10,:).^2))*1e5;
              
            % update temporary table with column info
            tmp_table.part_col=part_col'/du;
            tmp_table.tot_col=tot_col'/du;
            tmp_table.part_col_err=part_col_err'/du;
            tmp_table.tot_col_err=part_col_err'/du; % add error for sonde?? likely very small
            
            tmp_table.tot_col_smooth=tot_col_smooth'/du;
            tmp_table.tot_col_smooth_bk=tot_col_smooth_bk'/du;
            tmp_table.tot_col_smooth_pa=tot_col_smooth_pa'/du;
            
            tmp_table.part_col_smooth_bk=part_col_smooth_bk'/du;
            tmp_table.part_col_smooth_pa=part_col_smooth_pa'/du;
            
            tmp_table.col_below_rl=col_below_rl'/du;
            
            
        case 'NO2'

            % load NO2 data
            tg_num_dens=ncread(flist{i},'NO2_concentration'); % in mol/m3

            % convert to molec/cm3
            tg_num_dens=tg_num_dens*6.022140857e+17; 

% %             % get P, T
% %             pressure_hPa=ncread(flist{i},'pressure'); 
% %             temperature_K=ncread(flist{i},'temperature'); 
% %             % air number density (convert to molec/cm3)
% %             num_dens=((6.022e23*pressure_hPa.*100)./(8.314*temperature_K))*1e-6;           
% %             % get VMR
% %             vmr=NO2_num_dens./num_dens;
            
            
            % satellite partial column
            part_col=integrate(altitude*1e5, tg_num_dens(:,goodind),...
                               lowlim*1e5, highlim*1e5, 'midpoint');

            %%% AVK smoothing
            [lut_prof,lut_avk]=read_DOAS_prof_avk(2,[tmp_table.year,tmp_table.fractional_time],...
                                                  AVK_LUT_dir );

            % get avk smoothed partial column
            part_col_smooth=integrate_smooth(altitude*1e5,tg_num_dens(:,goodind),...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               [[0.5:59.5]*1e5]', lut_prof(:,3:end)', lut_avk(:,3:end)', 'interp');      
            % same for UV
            [lut_prof,lut_avk]=read_DOAS_prof_avk(3,[tmp_table.year,tmp_table.fractional_time],...
                                                  AVK_LUT_dir );

            % get avk smoothed partial column
            part_col_smooth_uv=integrate_smooth(altitude*1e5,tg_num_dens(:,goodind),...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               [[0.5:59.5]*1e5]', lut_prof(:,3:end)', lut_avk(:,3:end)', 'interp');      

            % bruker
            [bk_alt,bk_layer_h,bk_ap,bk_avk]=read_bruker_prof_avk(2,tmp_table.mjd2k);
            
            % get Bruker avk smoothed partial column 
            part_col_smooth_bk=integrate_smooth( altitude*1e5,tg_num_dens(:,goodind),...
                               lowlim*1e5,highlim*1e5,'midpoint',...
                               bk_alt'*1e5, bk_ap', bk_avk', 'interp',bk_layer_h*1e5);      

            % estimate error as 1e8 molec/cm3 in each layer (Chris Sioris)
            error=ones(size(tg_num_dens(:,goodind)))*1e8;
            part_col_err=sqrt(sum(error(lowlim-9:highlim-10,:).^2))*1e5;
            % part_col_err=NaN(size(part_col));
                           
            % update temporary table with column info
            tmp_table.part_col=part_col';
            tmp_table.tot_col=part_col';
            tmp_table.part_col_err=part_col_err';
            tmp_table.tot_col_err=part_col_err';
            
            tmp_table.tot_col_smooth=part_col_smooth';
            tmp_table.tot_col_smooth_uv=part_col_smooth_uv';
            tmp_table.tot_col_smooth_bk=part_col_smooth_bk';
            
            % find top (or bottom) altitude where data exists
            for j=1:length(goodind)
                
                tmp=~isnan(tg_num_dens(:,goodind(j)));
                tmp=altitude(tmp);
                top_alt=[top_alt, tmp(end)];
                
            end
            
    end

    % save coordinates, just in case
    tmp_table.lat=latitude(goodind);
    tmp_table.lon=longitude(goodind);
    
%     % add VMR profile
%     tmp_table.vmr=vmr(lowlim-9:highlim-10,goodind)';

    % add num_dens profile
    tmp_table.num_dens=tg_num_dens(lowlim-9:highlim-10,goodind)';

    % update final table
    osiris=[osiris; tmp_table];

end

% get rid of single precision
osiris.dist=double(osiris.dist);
osiris.lat=double(osiris.lat);
osiris.lon=double(osiris.lon);

if strcmp(tg,'NO2')
    % add scaling factors
    [scaleto_40,scaleto_42,scaleto_45]=pair_no2_lut(osiris,lowlim,highlim);

    osiris.scaleto_40=scaleto_40;
    osiris.scaleto_42=scaleto_42;
    osiris.scaleto_45=scaleto_45;
end
% mean(top_alt)

%% save results

% create new folder for radius around PEARL (500km is the default)
if comparison_range==500
    save_dir=data_dir;
else
    save_dir=[data_dir num2str(comparison_range) 'km/'];
    if ~exist(save_dir,'dir'), mkdir(save_dir), end
end

% save data
switch tg
    case 'O3'
        save([save_dir 'OSIRIS_v5p10_O3_table' comment '.mat'],'osiris')
    case 'NO2'
        save([save_dir 'OSIRIS_v6p00_NO2_table' comment '.mat'],'osiris')
end

cd(cur_dir);

% data=LoadHDFEOSFile( 'OSIRIS-Odin_L2-O3-Limb-MART_v5-07_2002m0311.he5' );
% data=LoadHDFEOSFile( 'OSIRIS-Odin_L2-NO2-Limb-Chalmers-DOAS-OE_v03-00_2004m0311.he5' );
end


function [scaleto_40,scaleto_42,scaleto_45]=pair_no2_lut(osiris,lowlim,highlim)


AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';    

% need to work in LUT's own directory
cur_dir=(pwd);
cd([AVK_LUT_dir 'no2_avk_lut_v2_0/']);

scaleto_40=NaN(length(osiris.lat),1);
scaleto_42=NaN(length(osiris.lat),1);
scaleto_45=NaN(length(osiris.lat),1);

% select wavelength
% from OSIRIS DOAS window, not sure for MART (relative values should be
% similar for other wavelengths
lambda=456; 
    
n=0;
for i=1:length(osiris.lat)
    
    % display progress info
    disp_str=['Generating scaling factor ', num2str(i), '/', num2str(length(osiris.lat))];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
    
    % write first input file
    fid = fopen('input_file_no2_avk.dat', 'w');
    fprintf(fid, '%s\n', '*Input file for NO2 AMF interpolation program');
    fprintf(fid, '%s\n', '*');
    fprintf(fid, '%s\n', '*Wavelength (350-550 nm) ?');
    fprintf(fid, '%s\n', num2str(lambda));
    fprintf(fid, '%s\n', '*Latitude (-90 (SH) to +90 (NH)) ?');
    fprintf(fid, '%s\n', num2str(osiris.lat(i)));
    fprintf(fid, '%s\n', '*Longitude (-180 (- for W) to +180 (+ for E)) ?');
    fprintf(fid, '%s\n', num2str(osiris.lon(i)));
    fprintf(fid, '%s\n', '*Ground albedo flag: 1 for Koelemeijer dscd_vecbase and 2 for albedo value defined by the user');
    fprintf(fid, '%s\n', '1');
    fprintf(fid, '%s\n', '*Ground albedo value (if albedo flag = 1, put -99)');
    fprintf(fid, '%s\n', '-99');
    fprintf(fid, '%s\n', '*Name of the file with SZA values for interpolation (less than 30 characters) ?');
    fprintf(fid, '%s\n', 'DAY_FILE.dat');
    fprintf(fid, '%s\n', '*Display flag (1: display the results on the screen; 0: dont display the results on the screen)');
    fprintf(fid, '%s\n', '0');
    fprintf(fid, '%s\n', '*Unit of the interpolated NO2 vertical profile: 0->VMR; 1->molec/cm3');
    fprintf(fid, '%s\n', '1');
    fclose(fid);

    % write second input file
    % year, fractional day (jan. 1, 00:00 = 0)
    fid = fopen('DAY_FILE.dat', 'w');
    fprintf(fid, '%.0f\t%.3f\n', ...
            [osiris.year(i), osiris.fractional_time(i)]);
    fclose(fid);

    % run LUT
    [status, result] = dos('wine no2_avk_interpolation_v2_0.exe', '-echo');

    %% read in avk and profile results
    % Format string for each line of text (both files):
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
    startrow=6;

    % no2 file
    fid = fopen('no2_prof_output.dat','r');
    % read data
    dataArray = textscan(fid, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,startrow-1, 'ReturnOnError', false);
    fclose(fid);
    % assign no2 profile variable
    lut_prof = [dataArray{1:end-1}];
    
    lut_prof(1:2)=[];
    
    % get scaling factors for various profile heights
    scaleto_40(i)=(sum(lut_prof(lowlim+1:40))*1e5) ./ ...
                  (sum(lut_prof(lowlim+1:highlim))*1e5);

    scaleto_42(i)=(sum(lut_prof(lowlim+1:42))*1e5) ./ ...
                  (sum(lut_prof(lowlim+1:highlim))*1e5);

    scaleto_45(i)=(sum(lut_prof(lowlim+1:45))*1e5) ./ ...
                  (sum(lut_prof(lowlim+1:highlim))*1e5);

              
%     % AVK file
%     fid = fopen('no2_avk_output.dat','r');
%     % read data
%     dataArray = textscan(fid, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,startrow-1, 'ReturnOnError', false);
%     fclose(fid);
%     % assign no2 avk variable
%     lut_avk = [dataArray{1:end-1}];

%     clearvars dataArray startrow formatSpec;    
    
end

cd(cur_dir)

fprintf('\n')
disp('Done')

end


