function reformat_OSIRIS_hdf( tg )
%REFORMAT_OSIRIS(tg): reformat OSIRIS HDF data as a matlab table


switch tg
    case {'O3','O3UV'}
        lowlim=14;
        highlim=52;
    case 'NO2'
        lowlim=14;
        highlim=40;
end


%% make list of data files

cur_dir=pwd;
data_dir='/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/';

cd([data_dir tg]);

tmp=dir('*.he5');
flist={tmp.name};
flist(1:2)=[];

%% setup output table
osiris=table;   


%% loop through files and read in/reformat data

year_prev=2000; % for progress display
du=2.687e16;

for i=1:length(flist)
    
    % load file (functions from OSIRIS website)
    try
        data=read_osiris_hdf(flist{i});
    catch
        disp(['Couldn''t read file ' flist{i}])
        continue
    end
    
    % get rid of annoying single precision
    data = structfun(@double, data, 'uniformoutput', 0);
    
    % filter by geolocation (only keep measurements within 500km of PEARL)
    range_km=dist_to_PEARL(data.Latitude,data.Longitude);
    goodind=find(range_km<=500);

    % skip if there are no measurements close to PEARL
    if isempty(goodind), continue, end
    
    % get time info
    % time is TAI 1993: seconds since 1993 Jan. 1, 00:00 (taken as equivalent to UTC)
    mjd2k=data.Time(goodind)/(24*3600)+ft_to_mjd2k(0,1993);
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
            
            % satellite partial column
            part_col=integrate(data.Altitude*1e5, data.O3NumberDensity(:,goodind),...
                               lowlim*1e5, highlim*1e5, 'trapez');
            
            % surface to 14km partial column from sonde data
            tot_col=NaN(size(part_col));

            for jj=1:length(part_col)

                [sonde_alt_grid,sonde_vmr,sonde_num_dens,~] =...
                    interp_ozonesonde( year(jj), fractional_time(jj) );

                sonde_conc=sonde_vmr.*sonde_num_dens;

                tmp=integrate(sonde_alt_grid*1e2,sonde_conc,2000,lowlim*1e5,'midpoint');
                tot_col(jj)=part_col(jj) + tmp;

            end
            
            % calculate error
            error=data.O3Precision(:,goodind).*data.RTModel_AirDensity...
                  (1:length(data.Altitude),goodind);
              
            % estimate column error (14-52 km limits)
            part_col_err=sqrt(sum(error(lowlim+1:highlim,:).^2))*1e5;
              
            % update temporary table with column info
            tmp_table.part_col=part_col'/du;
            tmp_table.tot_col=tot_col'/du;
            tmp_table.part_col_err=part_col_err'/du;
            tmp_table.tot_col_err=part_col_err'/du; % add error for sonde?? likely very small
            
        case 'NO2'
            

            % satellite partial column
            part_col=integrate(data.Altitude*1e5, data.NO2NumberDensity(:,goodind),...
                               lowlim*1e5, highlim*1e5, 'trapez');

            % calculate error
            error=data.NO2Precision(:,goodind).*data.RTModel_NeutralDensity...
                  (1:length(data.Altitude),goodind);
              
            % estimate error (14-40 km limits)
            part_col_err=sqrt(sum(error(lowlim/2+1:highlim/2+1,:).^2))*1e5;
              
            % update temporary table with column info
            tmp_table.part_col=part_col';
            tmp_table.tot_col=part_col';
            tmp_table.part_col_err=part_col_err';
            tmp_table.tot_col_err=part_col_err';
            
    end

    % save coordinates, just in case
    tmp_table.lat=data.Latitude(goodind);
    tmp_table.lon=data.Longitude(goodind);
    
    % update final table
    osiris=[osiris; tmp_table];

end

%% save results
switch tg
    case 'O3'
        save([data_dir 'OSIRIS_v5p07_O3_table.mat'],'osiris')
    case 'NO2'
        save([data_dir 'OSIRIS_v3p00_NO2_table.mat'],'osiris')
end

cd(cur_dir);

% data=LoadHDFEOSFile( 'OSIRIS-Odin_L2-O3-Limb-MART_v5-07_2002m0311.he5' );
% data=LoadHDFEOSFile( 'OSIRIS-Odin_L2-NO2-Limb-Chalmers-DOAS-OE_v03-00_2004m0311.he5' );
end

% % % %% functions from OSIRIS website
% % % function [data] = LoadHDFEOSFile( filename )
% % % 
% % %     swathinfo =  h5info(filename, '/HDFEOS/SWATHS');
% % %     variables = L_LoadGroupVariables( filename, swathinfo);
% % %     
% % %     nvar  = numel(variables);
% % %     data = [];
% % %     if (nvar > 0)
% % %         for i = 1:nvar
% % %             name = variables{i}.name;
% % %             data.(name) = variables{i}.value;
% % %             bad = (data.(name) == -9999.0);
% % %             data.(name)(bad) = NaN;
% % %         end
% % %     end
% % % end
% % % %%%
% % % function [variables] = L_LoadGroupVariables( filename, group)
% % % 
% % %     variables = {};
% % %     for i = 1: numel(group.Groups)
% % %         [v] = L_LoadGroupVariables( filename, group.Groups(i));   
% % %         variables = [variables; v];
% % %     end
% % %     
% % %     for i = 1: numel(group.Datasets)
% % %         v = L_LoadDataSet( filename, group.Name, group.Datasets(i) );
% % %         variables = [variables; v];
% % %     end
% % % end
% % % %%%
% % % function [v] = L_LoadDataSet( filename, groupname, dataset )
% % % 
% % %     entry       = [];
% % %     datasetname = [ groupname, '/',dataset.Name];
% % %     value       = h5read( filename, datasetname );
% % %     entry.name  = dataset.Name;
% % %     entry.value = value;
% % %     v           = {entry};
% % % end
% % % 
