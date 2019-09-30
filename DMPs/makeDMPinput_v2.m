function makeDMPinput_v2(varargin)
% This program will create the input for the DMP products created by Gloria
% Manney et al. 
%
% The user can input three variables: sitename, source_type, and instrument
% name. No inputs are necessary. Defaults are Eureka, NDACC, 125HR
%
% Sites must be one of the possiblesites listed below
% Available source types: NDACC HDF files, TCCON NetCDF files, .dat, AMES
%
% In addition, the 'general' option can be selected to calculate LoS
% parameters for a set of pre-defined SZA -- code finds available SZA for
% each day specified interval, and returns time, SZA and SAA
%
% Dan Weaver, University of Toronto
%
% Version 1.0: November 2015
% Version 2.0: August 2016, ensures that the latitute & longitude avoid
% values between -180 and -360. 
% "the program can accept -180 to 180 or 0 to 360" - Luis Millan, JPL
%
% modified by Kristof Bognar for use with GBS instruments, Nov. 2017
%
% Example: makeDMPinput_v2('Eureka','general','NO2_UV');

% addpath('DMPcode_v2');
%% OPTIONS
% choose what altitude grid to use
% if define_z = 0, it will use a 1-km grid up to 120 km
% if = 1, it will use the retrieval grid of the NDACC file
% if = 2, use a .mat file containing the altitude grid in variable "z".
% Define this file in the appropriate sub-program
define_z = 1;
% 
PutFilesInYearlyFolders = 1; % if ==1, will put files into yearly folders.
%% set up
possiblesites = {'Eureka', 'Thule','Kiruna','Pokerflat','Harestua','NyAlesund'};
possiblesources = {'bruker','TCCON','ames','dat','general','meas_time','osiris'};
instrument_name_override = 0;
%% OPTIONS ****************************************************************
switch nargin 
    case 0
        sitename = 'Eureka';
        source_type = 'bruker';
    case 1
        sitename = varargin{1};
        if ~strcmp(sitename,possiblesites) % check to ensure site aligns with expected sites
            error(['Input site ' sitename ' not recognized.'])
        end
        source_type = 'bruker';
    case 2
        %
        sitename = varargin{1};
        if ~strcmp(sitename,possiblesites) % check to ensure site aligns with expected sites
            error(['Input site ' sitename ' not recognized.'])
        end
        %
        source_type = varargin{2}; % NDACC, TCCON, dat, ames
        if ~strcmp(source_type,possiblesources)
            error(['Input source type ' source_type ' not recognized.'])
        end
    case 3
        %
        sitename = varargin{1};
        if ~strcmp(sitename,possiblesites) % check to ensure site aligns with expected sites
            error(['Input site ' sitename ' not recognized.'])
        end
        %
        source_type = varargin{2}; % NDACC, TCCON, dat, ames
        if ~strcmp(source_type,possiblesources)
            error(['Input source type ' source_type ' not recognized.'])
        end
        %
        instrument_name_override = 1; % force what instrument name is output in DMP input file
        tg_tag = varargin{3};
        %
end
% *************************************************************************
% outputdir = ['DMP_input_files/' sitename '/' source_type '/'];
outputdir = ['/home/kristof/work/DMP/DMP_input_files/' tg_tag '/'];
%
%% SET UP *****************************************************************
% 
Nmsmts = 0;
time = struct;
location = struct;
do_osiris=false;
%
%% IMPORT DATA
switch source_type
    case 'bruker'
        % NDACC only (quality controlled)
        % [time,location,sza,saa,z] = prepDMPinput_NDACC(sitename,define_z,tg_tag);
        % all measurements
        [time,location,sza,saa,z] = prepDMPinput_bruker_dat();
        
        LOSinfo = calc_LoS_bruker( z, sza, location, time );
        outputdir = ['/home/kristof/work/DMP/DMP_input_files/'];
    case 'TCCON'
        [time,location,sza,z] = prepDMPinput_TCCON(sitename);
    case 'ames'
        [time,location,sza,z] = prepDMPinput_NDACCames(sitename);
    case 'dat'
        [time,location,sza,z] = prepDMPinput_dat_v2(sitename);
    case 'general'
        % inputs for fixed list of sza every day (huge numbers)
        % do time/sza/saa generation and LoS calculation in one go
        [time,location,sza,saa,z] = sza_saa_for_LoS(1999, 2017, [72:2:90]);
        LOSinfo = calc_LoS_GBS( tg_tag, location, sza, saa, z, [72:2:90] );
    case 'meas_time'
        % inputs for actual GBS/SAOZ measurements
        % do time/sza/saa generation and LoS calculation in one go
        [time,location,sza,saa,z] = sza_saa_for_LoS_meas_time(tg_tag);
        LOSinfo = calc_LoS_GBS_interp( tg_tag, location, sza, saa, z );
    case 'osiris'
        % load osiris file for PEARL, reformat accordingly
        [ LOSinfo, time, z ] = calc_LoS_OSIRIS( tg_tag );
        do_osiris=false;
                        
end

% disp(['Read SZA and time info, ' num2str(length(sza)) ' datapoints found']);

%
%% Calculate the line-of-sight latitude and longitude *********************
% % % LOSinfo = calcLOSpts_v2(z,sza,location, time); % contains lat/lon info for input z points

% LOSinfo = calc_LoS_GBS( instr_name, location, sza, saa, z, [72:2:90] );

disp('Calculated LoS info for all datapoints');

%  
%% make DMP input file ****************************************************
% 
if instrument_name_override == 1
    switch source_type
        case 'bruker'
            instrument_name = ['125HR'];
        case 'general'
            instrument_name = ['GBS_' tg_tag];
        case 'meas_time'
            instrument_name = ['DOAS_' tg_tag];
        case 'osiris'
            instrument_name = ['OSIRIS_' tg_tag];
    end
    

else
%     [ instrument_name ] = NDACC_instrument( sitename );
end
%
Nmsmts = length(time.year);
Nz = length(z);
%

if ~do_osiris
    n=0;
    for j = 1:Nmsmts
        % organize output of files
        if PutFilesInYearlyFolders == 1 
            outputdirfinal = [outputdir num2str(time.year(j)) '/'];
            if ~exist(outputdirfinal,'dir')
                mkdir(outputdirfinal)
            end
        else
            outputdirfinal = outputdir;
            if ~exist(outputdirfinal,'dir')
                mkdir(outputdirfinal)
            end
        end

        % calculate the time of each measurement in seconds for filename
        timesecs = round(time.hour(j)*(60*60) + time.min(j)*60 + time.sec(j)); % round to avoid decimal that changes filename structure

        % filename
        strmonth = sprintf('%02d', time.month(j)); strday = sprintf('%02d', time.day(j));
        DMP_str=[sitename '_' instrument_name '_DMP_' num2str(time.year(j)) strmonth strday '.' num2str(timesecs) '.txt'];
        DMPfilename = [outputdirfinal DMP_str];

        % display progress info
        disp_str=['Writing file ' DMP_str];
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    


        if exist(DMPfilename,'file'), continue; end % avoid overwriting files
        % write DMP input file
        fileID = fopen(DMPfilename,'w');
        fprintf(fileID,'year day month time(sec) lat(deg) lon(deg) altitude(km) \r\n');
        for alt = 1:Nz
            fprintf(fileID,'%4d\t %02d\t %02d\t %5d\t %02.2f\t %03.2f\t %3.2f\t \r\n',...
                    time.year(j),time.day(j),time.month(j),timesecs,...
                    LOSinfo.Lat(alt,j),LOSinfo.Lon(alt,j),LOSinfo.z(alt));
        end
        fclose(fileID); 

    end
    
else % osiris yearly files
    
    outputdirfinal = outputdir;
    if ~exist(outputdirfinal,'dir')
        mkdir(outputdirfinal)
    end
    
    % calculate the time of each measurement in seconds for filename
    timesecs = round(time.hour*(60*60) + time.min*60 + time.sec); % round to avoid decimal that changes filename structure

    n=0;
    for yy=unique(time.year)'
    
        % find corresponding indices
        ind_yy=find(time.year==yy);
        
        % filename
        DMP_str=[sitename '_' instrument_name '_DMP_' num2str(yy) '.txt'];
        DMPfilename = [outputdirfinal DMP_str];
        
        % display progress info
        disp_str=['Writing file ' DMP_str];
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    
        
        if exist(DMPfilename,'file'), continue; end % avoid overwriting files
        % write DMP input file
        fileID = fopen(DMPfilename,'w');
        fprintf(fileID,'year day month time(sec) lat(deg) lon(deg) altitude(km) \r\n');
        for j=ind_yy'
            fprintf(fileID,'%4d\t %02d\t %02d\t %5d\t %02.2f\t %03.2f\t %3.2f\t \r\n',...
                    time.year(j),time.day(j),time.month(j),timesecs(j),...
                    LOSinfo.Lat(j),LOSinfo.Lon(j),LOSinfo.z);
        end
        fclose(fileID); 
        
        
    end
    
end
% good bye
fprintf('\n')
disp([num2str(Nmsmts) ' DMP files for ' sitename ' ' instrument_name ' have been created. Enjoy!']);
%

end
%
%%