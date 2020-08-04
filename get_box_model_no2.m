function get_box_model_no2()
%GET_BOX_MODEL_NO2(): model no2 columns for given years/days
%
% Input atmosphere is taken from ozonesonde measurements interpolated to
% the specified time
%   Sondes are extended to 60km using the US standard atmosphere
%
% NO2 columns are calculated for altitude limits specified in the code
%
% Surface albedo is calculated from the NDACC look-up tables (provided with
% the AMF LUTs)
%
% For details of the 'model_in' structure, see run_boxmodel.m
%
% Kristof Bognar, January 2018


% % % for highlim=[40,42,45]

%% input

% albedo options
% albedo_select='clim'; % climatology value for given month
albedo_select='interp'; % interpolated climatology (clim. at middle of month)
% albedo_select='fix'; % fixed at 0.9

% sonde options
% sonde_select='nearest'; % pick sonde profile nearest to given time
sonde_select=''; % interpolate sonde data to given time

% OSIRIS data starts in 2002, but there's no UT-GBS data that year
% run model for all the days
% startyear=2003;
startyear=1999;
endyear=2020;

% integration limits for NO2 column (satval paper)
% lowlim=12;
% highlim=40;
% max_days_in=300; % use February 15 - October 27, where max SZA is ~92.5

% integration limits for NO2 column (ozone depletion)
lowlim=10;
highlim=60;
max_days_in=151; % to end of May only

%% setup

% output table
model_no2=table;
model_no2.mjd2k=0;
model_no2.year=0;
model_no2.fractional_time=0;
model_no2.sza=0;
model_no2.no2=0;
model_no2.day=0;
model_no2.ampm=0;

% box model code directory
% addpath('/home/kristof/work/models/NO2_box_model/matlab')

% load standard atm
load('/home/kristof/work/NDACC/HDF4_data_submission/US_standard_AFGL1976/USstandard.mat');

% set constant model input
model_in.lat=80.05;
model_in.ndmr=1;
model_in.aero_scalefactor=1;
model_in.diurnalO3=0;

% define altitude grid
alt_grid=0.5:60.5;
grid_bound=[];
for iii=0:max(size(alt_grid))-1
    grid_bound=[grid_bound; [iii,iii+1]];
end

model_in.z=alt_grid';

% initialize input atmosphere array
model_in.atmos=[model_in.z zeros(length(model_in.z),6)];

%% loop through years
for i=startyear:endyear
    
    % pick dates and deal with leap years
    startday=46;
    if mod(i,4)==0, 
        max_days=max_days_in+1;
    else
        max_days=max_days_in;
    end
    
    disp(num2str(i));

    %% loop through days 
    n=0;
    nn=0;
    for j=startday:max_days

        % display progress info (with progress bar, cuz I was bored)
        percent=(j-(startday-1))*100/(max_days-(startday-1));
        disp_str=['Progress: [' repmat(' ',1,50) ']; day ' num2str(j)];
        disp_str(12:nn+11)='=';
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    
        
        if percent>2*(nn+1)
            nn=nn+1;
        end
        
        % set changing model input
        model_in.jday=j;
        
        switch albedo_select
            case 'clim' % albedo with fixed monthly values
            model_in.albedo=albedo_lookup(80.05, -86.4, i, j-1);
            
            case 'interp' % albedo interpolated, taking monthly values at middle of the month
            model_in.albedo=albedo_lookup(80.05, -86.4, i, j-1, 2);
            
            case 'fix' % constant albedo
            model_in.albedo=0.9;
        end
        
        %% get ozonesonde data
        % get data at noon of given day
        [ alt_grid_tmp, vmr_tmp, n_tmp, ~, ~, T_tmp]=interp_ozonesonde(i,j-0.5,sonde_select);
        
        alt_grid_tmp=alt_grid_tmp/1000;
        
        %% extend to 60 km using interpolated standard atm 
        
        % keep imported 10m grid, but extend to 60 km
        alt_grid_tmp_ext=[alt_grid_tmp;[40.015:0.01:60.005]'];
        
        % temporary arrays for individual profiles (extend NaNs to 60 km)
        T=NaN(size(alt_grid_tmp_ext));
        T2=NaN(size(alt_grid));
        T(1:4000)=T_tmp;
        
        vmr=NaN(size(alt_grid_tmp_ext));
        vmr2=NaN(size(alt_grid));
        vmr(1:4000)=vmr_tmp;
        
        num_dens=NaN(size(alt_grid_tmp_ext));
        num_dens2=NaN(size(alt_grid));
        num_dens(1:4000)=n_tmp;
        
        % find where profiles are missing
        ind_nan_T=find(isnan(T));
        ind_nan_vmr=find(isnan(vmr)); 
        ind_nan_num_dens=find(isnan(num_dens));
        
        % interpolate standard atm onto temporary altitude grid to
        % extend profile (units match)
        T(ind_nan_T)=interp1(alt_grid,USstandard_T,alt_grid_tmp_ext(ind_nan_T)); %K
        vmr(ind_nan_vmr)=interp1(alt_grid,USstandard_o3,alt_grid_tmp_ext(ind_nan_vmr));%unitless
        num_dens(ind_nan_num_dens)=interp1(alt_grid,USstandard_dens,alt_grid_tmp_ext(ind_nan_num_dens)); %cm-3

        % width of smooting window in km:
        fw_km=6; 
        % half width of window in grid points (10m grid)
        hw=fw_km*500/10; 
        
        % smooth out kinks
        T=boxcar_for_profiles(alt_grid_tmp_ext,T,hw,ind_nan_T(1)-1);
        vmr=boxcar_for_profiles(alt_grid_tmp_ext,vmr,hw,ind_nan_vmr(1)-1);
        num_dens=boxcar_for_profiles(alt_grid_tmp_ext,num_dens,hw,ind_nan_num_dens(1)-1);
        
        %% calculate effective layer values for model grid
        %(layers are too thick to just interpolate)
        for iii=1:max(size(alt_grid))
            
            % indices of layer boundaries
            ind_eff=find(alt_grid_tmp_ext>=grid_bound(iii,1) & alt_grid_tmp_ext<grid_bound(iii,2));
            
            % use normal mean, so if layer contains NaNs (i.e. profile doesn't cover 
            % entire layer) then result is NaN
            T2(iii)=nanmean(T(ind_eff));
            vmr2(iii)=nanmean(vmr(ind_eff));
            num_dens2(iii)=nanmean(num_dens(ind_eff));
            
        end
        
        % if sonde data is not extended with standard atm:
        % replace NaNs with zeros -- model recognizes 0 as no input
        % clarify model behavior when columns are partially filled
%         T2(isnan(T2))=0;
%         vmr2(isnan(vmr2))=0;
%         num_dens2(isnan(num_dens2))=0;
        
        % assign model input
        model_in.atmos(:,2:4)=[T2',num_dens2',[vmr2.*num_dens2]'];
        

        %% run box model and save output
        
        data=run_boxmodel(model_in);
        
        % remove second noon profile -- same as first
        data.LST(end)=[];
        data.SZA(end)=[];
        data.NO2(end,:)=[];
        
        % integrate NO2 data
        no2=integrate( alt_grid*1e5, data.NO2, lowlim*1e5, highlim*1e5, 'midpoint');
        
        % calculate fractional time in UTC
        % UTC from LST (LST is calculated from the position of the sun, the
        % model takes no longitude or timezone information)
        % leave numbers that are >24, so ft takes into account that utc
        % times slip to next day
        utc=data.LST+(86.4/15);

        % fractional time
        ft=data.jday-1+utc/24;
        
        % ampm flag
        ampm=zeros(size(data.LST));
        ampm(data.LST>=12)=1;
        
        % save results
        table_tmp=table;
        table_tmp.mjd2k=ft_to_mjd2k(ft,i);
        table_tmp.year=ones(size(ft))*i;
        table_tmp.fractional_time=ft;
        table_tmp.sza=data.SZA;
        table_tmp.no2=no2';
        table_tmp.day=ones(size(ft))*data.jday;
        table_tmp.ampm=ampm;
        
        % sort data by time
        [~,sortind]=sort(table_tmp.fractional_time);
        table_tmp=table_tmp(sortind,:);
        
        % add to final variable
        model_no2=[model_no2; table_tmp];

    end
    
    fprintf(repmat('\b',1,n+1));
    fprintf('\n')
    
end

% remove dummy row
model_no2(1,:)=[];

% save data
save(['/home/kristof/work/models/NO2_box_model/no2_cols_'...
      num2str(lowlim) '-' num2str(highlim) '_km.mat'], 'model_no2')

disp('Done')

% % % end

% pause(60)
% ! sudo shutdown -h now

end
