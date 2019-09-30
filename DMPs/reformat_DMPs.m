% function reformat_DMPs()
%REFORMAT_DMPS read DMPs from netcdf files and save as yearly matlab files


%% location of DMPs and saved files
% GBS fix SZA version
% dmp_dir='/home/kristof/aurora/eureka/ut-gbs/DMP/Eureka_GBS_DMPs_fixSZA/';
% save_dir='/home/kristof/work/DMP/DMP_GBS_fixSZA/';
% species='NO2_VIS';

% GBS and SAOZ measurement times
% dmp_dir='/home/kristof/aurora/eureka/ut-gbs/DMP/Eureka_BRUKER_DOAS_OSIRIS_DMPs/';
% save_dir='/home/kristof/work/DMP/DMP_DOAS_meas_time/';
% species='DOAS_O3_VIS';

% bruker 2006-2015 (Dan's files on deluge)
dmp_dir='/home/kristof/deluge/projects/PEARL_FTS/DMPs/Eureka_NDACC/';
save_dir='/home/kristof/work/DMP/DMP_bruker/';
species='';

% bruker 2016-2017 
% dmp_dir='/home/kristof/aurora/eureka/ut-gbs/DMP/Eureka_BRUKER_DOAS_OSIRIS_DMPs/';
% save_dir='/home/kristof/work/DMP/DMP_bruker/';
% species='BRUKER';

% OSIRIS
% dmp_dir='/home/kristof/aurora/eureka/ut-gbs/DMP/Eureka_BRUKER_DOAS_OSIRIS_DMPs/';
% save_dir='/home/kristof/work/DMP/DMP_OSIRIS/';
% species='OSIRIS_O3';

% PARIS
dmp_dir='/home/kristof/work/DMP/DMP_PARIS/';
save_dir='/home/kristof/work/DMP/DMP_PARIS/';
species='';

cur_dir=pwd();

%% loop over years
for year=1999:2017
% for year=2006:2006
   
    % get list of files
    file_dir=[dmp_dir species '/' num2str(year) '/'];
    
    if ~exist(file_dir,'dir')
        continue
    else
        cd(file_dir)
    end
    
    flist=get_file_list(file_dir,'nc4');
    
    % initalize variables
    dmp_all=[];
    fractional_time=[];
    
    % check for empty files
    ind_empty=[];
    for i=1:length(flist)
        tmp=dir(flist{i});
        if tmp.bytes==0
            ind_empty=[ind_empty,i];
        end
    end
    if ~isempty(ind_empty), flist(ind_empty)=[]; end
    
    
    %% loop over files
    n=0;
    for i=1:length(flist)
        
        % display progress info
        disp_str=['Reading file ', num2str(i), '/', num2str(length(flist))];
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    
        
        % read date info (yyyymmdd, in a cell)
        % ncread doesn't reconize data type; need to read as HDF
        try
            date_cell=h5read(flist{i},'/dynprofs/Date');
            yyyy=date_cell{1}(1:4);
            mm=date_cell{1}(5:6);
            dd=date_cell{1}(7:8);
        catch
            % sometimes h5read fails, use file name then
            tmp=strsplit(char(flist{i}),'.');
            yyyy=tmp{1}(end-7:end-4);
            mm=tmp{1}(end-3:end-2);
            dd=tmp{1}(end-1:end);
        end
        
        % read time info (fractional hour)
        fh=ncread(flist{i},'dynprofs/Time');
        
        if length(unique(fh))~=1, error([flist{i} 'contains multiple times']); end
        
        HH=floor(fh(1));
        MM=floor((fh(1)-HH)*60);
        SS=round((fh(1)-HH-MM/60)*3600);
        
        % convert to fractional time and save
        [ft]=fracdate([yyyy '/' mm '/' dd ' ' num2str(HH) ':' num2str(MM) ':' num2str(SS)],...
                        'yyyy/mm/dd HH:MM:SS');
                    
        fractional_time=[fractional_time, ft];

        
        % load data into a table
        dmp=table;
        dmp.lat=ncread(flist{i},'dynprofs/Lat');
        dmp.lon=ncread(flist{i},'dynprofs/Lon');
        dmp.alt=ncread(flist{i},'dynprofs/Altitude')';
        dmp.pres=ncread(flist{i},'dynprofs/Pressure')';
        dmp.temperature=ncread(flist{i},'dynprofs/Temperature')';
        dmp.theta=ncread(flist{i},'dynprofs/Theta')';
        dmp.pv=ncread(flist{i},'dynprofs/PV')';
        dmp.spv=ncread(flist{i},'dynprofs/sPV')';
        
        % save table in cell array
        dmp_all{i}=dmp;
        
        
    end
    
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    
%     save([save_dir species '_DMP_table_' num2str(year)], 'dmp_all','fractional_time')
    save([save_dir 'PARIS_DMP_table_' num2str(year)], 'dmp_all','fractional_time')
    
    fprintf('\n');
    fprintf(['Done with ' num2str(year) '\n']);
    
    
end

cd(cur_dir)

% end

