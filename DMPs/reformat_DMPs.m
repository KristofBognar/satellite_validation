function reformat_DMPs(type_in, local_data, loc_in, save_dir, new_only)
%REFORMAT_DMPS read DMPs from netcdf files and save as yearly matlab files
%
% Read DMP data from netCDF files (one per measurement) and save in matlab
% format
%
% INPUT:
%   type_in: instrument or dataset, see instr_list below
%   local_data:
%       true: data on local system
%       false: data on aurora (might be extremely slow, each variable requires
%              a new read command)
%   loc_in: 
%       if local data: path to DMP folder
%       if ~local data: local_datapath to aurora folder on local system (must be full path,
%                       or exist might mess up)
%   save_dir: where to save the data (must be full path)
%   new_only (optional):
%       true: skip years that have been saved before (default)
%       false: process and save all available years
%
%@Kristof Bognar
%   original version (2018): script to manually select what instrument to process
%   major update (2020): function to read any DMP based on input

%% input

% check instrument name
instr_list={'DOAS_O3_VIS','DOAS_NO2_VIS','DOAS_NO2_UV',...
            'BRUKER',...
            'BREWER69','BREWER192','BREWER223',...
            'OSIRIS_O3','OSIRIS_NO2'};
        
if ~any(strcmp(type_in,instr_list)) 
    error(['Instrument ' type_in ' not recognized'])
end

% redo all files, or skip files that already exist?
if nargin==4, new_only=1; end

%% location of DMPs

if local_data
    
    aurora_dir='';
    loc_in=[loc_in type_in '/'];
    
else
    
    % DOAS data
    if strfind(type_in,'DOAS')
    %     aurora_dir=['ground/eureka/gbs/dmp/' type_in '/'];
        aurora_dir=[type_in '/'];
    end

    if strfind(type_in,'BRUKER')
        % Dan's files (up to 2015)
        % 2018 processing
        % 2020 processing
    end

    if strfind(type_in,'BREWER')
        aurora_dir=['ground/eureka/brewer/dmp/' type_in '/'];
    end

    if strfind(type_in,'OSIRIS')
        aurora_dir=['ground/eureka/gbs/dmp/tmp_archive_all/Eureka_2018/' type_in '/'];
    end

    % % % bruker 2006-2015 (Dan's files on deluge)
    % % dmp_dir='/home/kristof/deluge/projects/PEARL_FTS/DMPs/Eureka_NDACC/';
    % % save_dir='/home/kristof/work/DMP/DMP_bruker/';
    % % species='';
    % % 
    % % % bruker 2016-2017 
    % % % dmp_dir='/home/kristof/aurora/eureka/ut-gbs/DMP/Eureka_BRUKER_DOAS_OSIRIS_DMPs/';
    % % % save_dir='/home/kristof/work/DMP/DMP_bruker/';
    % % % species='BRUKER';
    % % 
    % % cur_dir=pwd();
end

%% loop over all possible years
for yr=1999:year(now)
   
    % DMP file version
    if yr<2018 && ~strfind(type_in,'BREWER')
        file_ver=2018; % 2018 DMP processing; one file per meas.
    else
        file_ver=2020; % 2020 DMP processing; 3 files per meas.
    end
    
    % saved file
    save_name=[type_in '_DMP_table_' num2str(yr) '.mat'];
    
    % get list of files
    file_dir=[loc_in aurora_dir num2str(yr) '/'];
    
    % check if year exists in DMP output
    if ~exist(file_dir,'dir')
        % no: write warning
        disp(['No DMPs for ' num2str(yr) ' ' type_in ' data'])
        continue
    else
        % yes: check if we want to process this year
        if (new_only && exist([save_dir save_name],'file'))
            % don't overwrite old files
            disp(['DMPs already processed for ' num2str(yr) ' ' type_in ' data'])
            continue
        end
    end
    
    % list if DMP fles for given year
    if file_ver==2018 % all variables in one file
        flist=get_file_list(file_dir,'nc4');
    elseif file_ver==2020 % variable groups broken up, need dyneql file
        flist=get_file_list(file_dir,'nc4','*DynEqL*');
    end
    
    % initalize variables
    dmp_all=[];
    fractional_time=[];
    
    % check for empty files
    ind_empty=[];
    for i=1:length(flist)
        tmp=dir([file_dir flist{i}]);
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
        
        % read DMP file
        [ft,dmp]=read_DMP_files([file_dir flist{i}],file_ver);
        
        % save results
        fractional_time=[fractional_time, ft];
        dmp_all{i}=dmp;
        
    end
    
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    
    save([save_dir save_name], 'dmp_all','fractional_time')
    
    fprintf('\n');
    fprintf(['Done with ' num2str(yr) '\n']);
    
    
end

end

%%
function [ft,dmp]=read_DMP_files(fname,file_ver)

    if file_ver==2018
        
        % DMP files from 2018 processing (and before): single file for all data
        var_grp='dynprofs/';
        
    elseif file_ver==2020
        
        % DMP files from 2020 processing: separate dyn, jet, trop files
        var_grp='';
        
    end

    % read date info (yyyymmdd, in a cell)
    % ncread doesn't reconize strings; need to read as HDF
    
    try
        date_cell=h5read(fname,['/' var_grp 'Date']);
        yyyy=date_cell{1}(1:4);
        mm=date_cell{1}(5:6);
        dd=date_cell{1}(7:8);
    catch
        % sometimes h5read fails, use file name then
        tmp=strsplit(char(fname),'.');
        yyyy=tmp{1}(end-7:end-4);
        mm=tmp{1}(end-3:end-2);
        dd=tmp{1}(end-1:end);
    end
    
    % read time info (fractional hour)
    if file_ver==2018
        
        % old files: stored as number
        fh=ncread(fname,'dynprofs/Time');
        if length(unique(fh))~=1
            error([fname 'contains multiple times']);
        else
            fh=fh(1);
        end
        
    elseif file_ver==2020
        
        % new files: stored as string, same issue as date
        try
            fh=h5read(fname,'/Hour');
            if str2double(fh{1})~=str2double(fh{end})
                error([fname 'contains multiple times']);
            else
                fh=str2double(fh{1});
            end
        catch
            % sometimes h5read fails, use file name then
            tmp=strsplit(char(fname),'.');
            tmp=strsplit(char(tmp{2}),'_');
            fh=str2num(tmp{1})/3600;
        end
        
    end

    HH=floor(fh);
    MM=floor((fh-HH)*60);
    SS=round((fh-HH-MM/60)*3600);

    % convert to fractional time and save
    [ft]=fracdate([yyyy '/' mm '/' dd ' ' num2str(HH) ':' num2str(MM) ':' num2str(SS)],...
                    'yyyy/mm/dd HH:MM:SS');


    % load data into a table
    
    dmp=table;
    dmp.lat=ncread(fname,[var_grp 'Lat']);
    dmp.lon=ncread(fname,[var_grp 'Lon']);
    dmp.alt=ncread(fname,[var_grp 'Altitude'])';
    dmp.pres=ncread(fname,[var_grp 'Pressure'])';
    dmp.temperature=ncread(fname,[var_grp 'Temperature'])';
    dmp.theta=ncread(fname,[var_grp 'Theta'])';
    dmp.pv=ncread(fname,[var_grp 'PV'])';
    dmp.spv=ncread(fname,[var_grp 'sPV'])';
    
end

