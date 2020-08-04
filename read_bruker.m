function [] = read_bruker(tg,instrument)
%read_bruker() read bruker data from NDACC HDF files and save continuous timeseries
%   Save coth total columns and

%% setup
satval=0;
if nargin==1, instrument='bruker'; end % can be set to paris

% check tg input
valid_tg={'O3','NO2','HCl','HNO3','ClONO2','HF'};
if ~any(strcmp(valid_tg,tg)), error('Select valid trace gas'); end


% location of bruker files
if satval
    filedir=['/home/kristof/work/' instrument '/' tg '/'];
    save_path=['/home/kristof/work/' instrument '/'];
else

    if ~strcmpi(instrument,'bruker'), error('Set file path'), end
    
    server_path='/home/kristof/atmosp_servers/';
    save_path='/home/kristof/work/bruker/PEARL_ozone_depletion/';
    
    switch tg
        case 'O3'
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDF_CAMS/'];
        case 'NO2'
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDF_NO2/no2.error.fixed/'];
        otherwise
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDFs_for_plotting/' lower(tg)];
    end
    
end

cur_dir=pwd();

cd(filedir)
temp = dir(['groundbased*' lower(tg) '*.hdf']); 
flist = {temp.name}; % cell array of file names

% variables to store data
time_mjd2k=[];
column=[];
err_sys=[];
err_rand=[];
sza=[];
saa=[];

failed={};

%% loop over HDF files and read data
n=0;
for i=1:length(flist)

    fname=flist{i};

    % display progress info
    disp_str=['Reading ' fname];
    % stuff to delete last line and reprint updated message
    fprintf(repmat('\b',1,n));
    fprintf(disp_str);
    n=numel(disp_str);    
        
    % read variables
    try
        time_mjd2k=[time_mjd2k, hdfread(fname,'DATETIME')];
        column=[column, hdfread(fname,[tg '.COLUMN_ABSORPTION.SOLAR'])];
        err_rand=[err_rand, hdfread(fname,...
            [tg '.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.RANDOM.STANDARD'])];
        err_sys=[err_sys, hdfread(fname,...
            [tg '.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.SYSTEMATIC.STANDARD'])];
        sza=[sza, hdfread(fname,'ANGLE.SOLAR_ZENITH.ASTRONOMICAL')];
        saa=[saa, hdfread(fname,'ANGLE.SOLAR_AZIMUTH')];
    catch
        failed{end+1}=fname;
    end
end

fprintf('\n')

if ~isempty(failed)
    disp('Could not read the following files:')
    for i=1:length(failed)
        disp(failed{i})
    end
end

% convert to date
dates=mjd2k_to_date(time_mjd2k');

[ft,years]=fracdate(dates);

months=month(dates);
days=day(dates);

% get rid of annoying single precision
column=double(column);
err_sys=double(err_sys);
err_rand=double(err_rand);
sza=double(sza);
saa=double(saa);

% convert to table
bruker=table(time_mjd2k',years,months,days,ft,column',err_sys',err_rand',sza',saa',dates);

bruker.Properties.VariableNames={'mjd2k', 'year', 'month', 'day', 'fractional_time',...
                                 'tot_col', 'tot_col_err_sys', 'tot_col_err_rand',...
                                 'sza', 'saa' , 'DateTime'};

bruker.Properties.VariableUnits={'','','','','jan 1, 00:00 = 0',...
                                 'molec/cm2','molec/cm2','molec/cm2','','',''};                             

                      
if satval

    %% match Bruker NO2 to initial submission in order to use the exact same 
    %% dataset in the revised paper
    if strcmp(tg,'NO2')

        disp('Matching Bruker NO2 to dataset in submitted version of satval paper')

        load(['/home/kristof/work/bruker/bruker_no2_submission_initial.mat'])

        [~,tmp]=setdiff(bruker.mjd2k,subm_init.mjd2k);

        bruker(tmp,:)=[];
    end

    %% partial columns/profiles

    num_dens=[];
    num_dens_native=[];
 
    switch tg
        case 'O3'

            % calculate 14-52 km partial columns for satellite comparisons
            [alt_bk,layer_h,~,~,part_prof,dof]=...
                read_bruker_prof_avk(tg,save_path,bruker.mjd2k,instrument);
            part_col_tmp=integrate_nonuniform(...
                alt_bk*1e5,part_prof,14e5,52e5,'midpoint', layer_h*1e5 );

            % save number density profiles
            for i=1:length(bruker.mjd2k)
                % interpolate to sat grid
                num_dens(i,:)=interp1(alt_bk,part_prof(i,:),14.5:51.5);
            end

            num_dens_native=part_prof(:,20:39);

            % add partial column to data
            bruker.part_col=part_col_tmp';

        case 'NO2'

            % calculate 12-40 km partial columns for satellite comparisons
            [alt_bk,layer_h,~,~,part_prof,dof,avk_trace]=...
                read_bruker_prof_avk(tg,save_path,bruker.mjd2k,instrument);
            part_col_tmp=integrate_nonuniform(...
                alt_bk*1e5,part_prof,12e5,40e5,'midpoint', layer_h*1e5 );

            % also save 12-60 km partial columns, matching DOAS range
            part_col_tmp_doas=integrate_nonuniform(...
                alt_bk*1e5,part_prof,12e5,60e5,'midpoint', layer_h*1e5 );

            % save number density profiles
            for i=1:length(bruker.mjd2k)
                % interpolate to sat grid
                num_dens(i,:)=interp1(alt_bk,part_prof(i,:),12.5:39.5);
            end
            
            % add partial column to data
            bruker.part_col=part_col_tmp';
            bruker.part_col_doas=part_col_tmp_doas';

            % save 12-40 km partial column dofs (12.205-40.17 km using bruker altitude range)
            bruker.dof_part=sum(avk_trace(:,17:36),2);

    end

    % save degrees of freedom (total column)
    bruker.dof=dof;

    % add num dens profile
    bruker.num_dens=num_dens;
    if strcmp(tg, 'O3'), bruker.num_dens_native=num_dens_native; end
    %%% no2 fine to here                             

    %% smooth using DOAS AVKs and save
    switch tg
        case 'O3'
            bruker=smooth_bruker( tg, save_path, bruker, instrument );
            if strcmp(instrument,'PARIS')
                paris=bruker;
                save([save_path instrument '_o3.mat'], 'paris');
            else
                save([save_path instrument '_o3.mat'], 'bruker');
            end
        case 'NO2'
            bruker=smooth_bruker( tg, save_path, bruker, instrument );
            save([save_path instrument '_no2.mat'], 'bruker');
    end
    
else
    save([save_path 'bruker_' lower(tg) '.mat'], 'bruker');
end


cd(cur_dir)

end

