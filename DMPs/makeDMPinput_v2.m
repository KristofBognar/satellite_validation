function makeDMPinput_v2(varargin)
% This program will create the input for the DMP products created by Gloria
% Manney et al. 
%
% The user can input four variables: sitename, source_type, trace gas tag,
% and the required years
%
% Sites must be one of the possiblesites listed below
%
% source types: 
%   original from Dan: NDACC HDF files, TCCON NetCDF files, .dat, AMES
%   'bruker': new wrapper for dat input
%   GBS options:
%       'general': fixed set of SZA, don't use
%       'meas_time': one DMP input for each measurements
%   'osiris': vertical DMP profile at OSIRIS tangent point
%
% param_in: GBS only: 'O3_VIS', 'NO2_VIS', 'NO2_UV'
%           brewer only: brewer instrument number, as string
%
% year input (string), works for 'bruker', 'meas_time', 'brewer': 
%   'all' to process all data
%   '2018' to process given year
%   '2018+' to process given year and all subsequent years
%   '2018-2020' to process a range of years (limits included)
%
% Dan Weaver, University of Toronto
% Version 1.0: November 2015
% Version 2.0: August 2016, ensures that the latitute & longitude avoid
% values between -180 and -360. 
% "the program can accept -180 to 180 or 0 to 360" - Luis Millan, JPL
%
% modified by Kristof Bognar:
% DMPs for the GBSs, Nov. 2017
% added support for yearly DMPs, May 2020
% DMPs for Brewers, May 2020
%
% 
%
% Example: makeDMPinput_v2('Eureka','meas_time','NO2_UV');

%% User guide

% Bruker data: use most up-to-date spDB_eur_<year>.dat file from Bruker team
%       convert the text file to matlab format and save
%       modify prepDMPinput_bruker_dat.m to so it loads the saved file
%       !! the saved files must contain the ALTITUDE variable (the
%           retrieval grid) from the Bruker HDF files (low to high elev in km)

% GBS data: use meas_time option (DMP inputs for each measurement time/SZA)!
%       create *_unfiltered files for each year using
%           merge_reanalysis_VCDs.m -- includes failed VCDs as well, in
%           case retrieval changes later
%       modify sza_saa_for_LoS_meas_time.m so it points to new *_unfiltered files
%       !! the standard NDACC altitude grid is used



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
possiblesources = {'bruker','TCCON','ames','dat','general','meas_time','osiris','brewer'};
possiblesources_nargin2 = {'bruker','TCCON','ames','dat'};
instrument_name_override = 0;
%% OPTIONS ****************************************************************
switch nargin 
%     case 0
%         sitename = 'Eureka';
%         source_type = 'bruker';
%     case 1
%         sitename = varargin{1};
%         if ~strcmp(sitename,possiblesites) % check to ensure site aligns with expected sites
%             error(['Input site ' sitename ' not recognized.'])
%         end
%         source_type = 'bruker';
    case 2
        %
        sitename = varargin{1};
        if ~strcmp(sitename,possiblesites) % check to ensure site aligns with expected sites
            error(['Input site ' sitename ' not recognized.'])
        end
        %
        source_type = varargin{2}; % NDACC, TCCON, dat, ames
        if ~strcmp(source_type,possiblesources_nargin2)
            error(['Input source type ' source_type ' not recognized.'])
        end
        %
        yr_out='all';
        %
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
        param_in = varargin{3};
        yr_out='all';
    case 4
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
        param_in = varargin{3};
        %        
        yr_out=varargin{4};
        %
    otherwise
        error('Need sitename, source_type, and tg_tag (see documentation)')
end
% *************************************************************************
outputdir = ['DMP_input_files/' sitename '/' source_type '/'];
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
    case 'bruker' % added by Kristof
        % NDACC only (quality controlled)
        % [time,location,sza,saa,z] = prepDMPinput_NDACC(sitename,define_z,tg_tag);
        
        % all measurements
        [time,location,sza,saa,z] = prepDMPinput_bruker_dat(yr_out);
        
        LOSinfo = calc_LoS_bruker( z, sza, location, time );
        outputdir = ['/home/kristof/work/DMP/DMP_input_files/BRUKER/'];

    case 'brewer' % added by Kristof
        % check input
        if isnan(str2double(param_in)), error('Provide Brewer number as string'), end
            
        % all measurements
        [time,location,sza,z] = prepDMPinput_brewer(str2double(param_in), yr_out);
        
        LOSinfo = calc_LoS_bruker( z, sza, location, time ); % use bruker code, same DS meas.
        outputdir = ['/home/kristof/work/DMP/DMP_input_files/BREWER' param_in '/'];
        
    case 'TCCON' % Dan's version
        [time,location,sza,z] = prepDMPinput_TCCON(sitename);
    case 'ames' % Dan's version
        [time,location,sza,z] = prepDMPinput_NDACCames(sitename);
    case 'dat' % Dan's version
        [time,location,sza,z] = prepDMPinput_dat_v2(sitename);
    case 'general' % added by Kristof
        % inputs for fixed list of sza every day (huge numbers)
        % do time/sza/saa generation and LoS calculation in one go
        error('Do not use')
        [time,location,sza,saa,z] = sza_saa_for_LoS(1999, 2017, [72:2:90]);
        LOSinfo = calc_LoS_GBS( param_in, location, sza, saa, z, [72:2:90] );
        outputdir = ['/home/kristof/work/DMP/DMP_input_files/GBS_' param_in '/'];
        
    case 'meas_time' % added by Kristof
        % inputs for actual GBS/SAOZ measurements
        % do time/sza/saa generation and LoS calculation in one go
        [time,location,sza,saa,z] = sza_saa_for_LoS_meas_time(param_in,yr_out);
        
        LOSinfo = calc_LoS_GBS_interp( param_in, location, sza, saa, z );
        
        outputdir = ['/home/kristof/work/DMP/DMP_input_files/DOAS_' param_in '/'];
        
    case 'osiris' % added by Kristof
        % load osiris file for PEARL, reformat accordingly
        [ LOSinfo, time, z ] = calc_LoS_OSIRIS( param_in );
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
        case 'brewer'
            instrument_name = ['BREWER' param_in];
        case 'general'
            instrument_name = ['GBS_' param_in];
        case 'meas_time'
            instrument_name = ['DOAS_' param_in];
        case 'osiris'
            instrument_name = ['OSIRIS_' param_in];
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