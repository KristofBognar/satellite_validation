% function reformat_SAOZ()
%REFORMAT_SAOZ reformat saoz data to match GBS file frmat
%   Warning: no time information exits in SAOZ files, only day and am/pm

%%% Use this to convert imported array to a table
%%% T = array2table(O3NO220052017,'VariableNames',{'Year','Month','Day','DoY','O3sr','O3ss','dO3sr','dO3ss','NO2sr','NO2ss','dNO2sr','dNO2ss'});


% load data (already saved as a table)
% % load('/home/kristof/work/SAOZ/O3_NO2_2005-2017.mat')
load('/home/kristof/work/SAOZ/O3_NO2_2005-2020.mat')

% convert to array for easier rearrangement
tmp=table2array(o3_no2);

% read dSCDs to add sza and calculate more precise time for twilight columns
% true: time and SZA are calculated from dscds (highest available 5deg SZA range)
%       used for Bognar et al. 2019
% false: time is from AMF file, and SZA is calculated using matlab
%        AMF tme is usually for SZA=90
read_dscd_time=false;

% reorganize by trace gas and am/pm
for tg=1:2
    
% %     % Bognar et al. 2019    
% %     if tg==1
% %         load('/home/kristof/work/SAOZ/AMF_O3_2005-2017.mat')
% %     else
% %         load('/home/kristof/work/SAOZ/AMF_NO2_2005-2017.mat')
% %     end
    if tg==1
        load('/home/kristof/work/SAOZ/AMF_O3_2005-2020.mat')
    else
        load('/home/kristof/work/SAOZ/AMF_NO2_2005-2020.mat')
    end
    
    % split sunrise and sunset columns into separate rows; add column for
    % ampm flag, add extra column of zeros for error calc (to mirror GBS cormat)
    % also add hour of the day from amf files (same number of rows)
    saoz=[ [tmp(:,1:3:4)  zeros(size(tmp,1),1)...
            tmp(:,4*tg+1:2:4*tg+3)  zeros(size(tmp,1),1)  amf_time(:,5) ];...
           [tmp(:,1:3:4)  ones(size(tmp,1),1)...
            tmp(:,4*tg+2:2:4*tg+4)  zeros(size(tmp,1),1)  amf_time(:,6) ]...
          ];

    % sort by year, doy, ampm
    saoz=sortrows(saoz,[1,2,3]);
    
    % calculate fractional time and julian date
    saoz(:,8)=saoz(:,2)-1+saoz(:,7)/24;
    saoz(saoz(:,3)==1 & saoz(:,7)<6,8)=saoz(saoz(:,3)==1 & saoz(:,7)<6,8)+1;    
    saoz(:,9)=ft_to_mjd2k(saoz(:,8),saoz(:,1));
    
    % remove temporary hour column
    saoz(:,7)=[];
    
    % convert to table
    saoz=array2table(saoz);
    saoz.Properties.VariableNames={'year' 'day' 'ampm'...
                                   'mean_vcd' 'std_vcd' 'sigma_mean_vcd'...
                                   'fractional_time' 'mjd2k'};
    % remove NaNs
    saoz(isnan(saoz.mean_vcd),:)=[];
    saoz(isnan(saoz.mjd2k),:)=[];
    
    %% add SZA/SAA info
    
    dscd_time=NaN(size(saoz.year));
    dscd_sza=NaN(size(saoz.year));

    if read_dscd_time % time and SZA from dSCD files
        
        disp('Reading time and SZA from dSCD file')
        load('/home/kristof/work/SAOZ/SAOZ_dSCD.mat')

        for i=1:length(saoz.year)

            % find matching twilight
            ind=find(dscd.year==saoz.year(i) & dscd.day==saoz.day(i) & dscd.ampm==saoz.ampm(i));

            if isempty(ind), continue, end

            % use max available SZA if below 91
            max_sza=min(max(dscd.sza(ind)),91);

            % find 5deg window
            ind2=find(dscd.sza(ind)<=max_sza & dscd.sza(ind)>=max_sza-5);

            % indices to use in time average
            ind_time=ind(ind2);

            % mean time
            dscd_time(i)=mean(dscd.fractional_time(ind_time));
            dscd_sza(i)=mean(dscd.sza(ind_time));

        end
        
    else
        disp('No dSCD file; SZA/SAA are calculated from AMF time')
    end
    
    % fill in time where dscds are missing
    % all dscds are missing if read_dscd_time=false
    ind_tmp=find(isnan(dscd_time));
    if ~isempty(ind_tmp)
        dscd_time(ind_tmp)=saoz.fractional_time(ind_tmp);
    else
        disp('All times found in dSCD file')
    end

    saoz.fractional_time=dscd_time;
    saoz.mjd2k=ft_to_mjd2k(saoz.fractional_time,saoz.year);

    % add SZA
    saoz.sza=dscd_sza;

    % Add SAA and fill in SZA where dscds are missing
    date_tmp=mjd2k_to_date(saoz.mjd2k);
    [az_tmp, sza_tmp] = SolarAzEl(date_tmp,80.053*ones(size(saoz.mjd2k)),...
                                  -86.416*ones(size(saoz.mjd2k)),...
                                  0.6*ones(size(saoz.mjd2k)));
    sza_tmp=90-sza_tmp;

    saoz.sza(ind_tmp)=sza_tmp(ind_tmp);
    saoz.saa=az_tmp;

    %% save
    
    if tg==1 % ozone
        save('/home/kristof/work/SAOZ/saoz_o3.mat','saoz');
    else
        saoz.mean_vcd=saoz.mean_vcd*1e15;
        saoz.std_vcd=saoz.std_vcd*1e15;
        save('/home/kristof/work/SAOZ/saoz_no2.mat','saoz');
    end
    
    
end


% end

