function [time,location,SZA,SAA,z] = prepDMPinput_NDACC(sitename,define_z,tag)
% This program will prepare the input for the DMP products from standard NDACC archived files
%
% November 2015
% Dan Weaver, University of Toronto

% % if exist('M:','dir') == 7
% %     data_drive = 'M:';
% % elseif exist('X:','dir') == 7
% %     data_drive = 'X:';
% % else
% %     error('no data drives found!')
% % end
% % %
% % addpath([data_drive 'Data/NDACC/'])  % code relies on the NDACC archive's HDF files
%
% Modified by Kristof Bognar to read GBS data from HDF format

% hdf_path='/home/kristof/work/DMP/HDF_tmp/';

%% OPTIONS
% choose what altitude grid to use
% if define_z = 0, it will use a 1-km grid up to 120 km
% if = 1, it will use the retrieval grid of the NDACC file
% if = 2, it will use a .mat file which must simply contain the altitude grid in a variable "z".
ALTfile = 'grid61.mat';
%
if define_z == 0
    z = 1:1:100;
elseif define_z == 2
    load(ALTfile);
end
%% SET UP *****************************************************************
%
N = 0;
Nmsmts = 0;
time = struct;
location = struct;
datastart = 0;
dataend = 0;

%
%% IMPORT DATA
%
% create list of all relevant files
% % filelist = ls([data_drive 'Data/NDACC/' sitename '/*.hdf']);
% % [Nfiles,~] = size(filelist);

switch tag
    case 'NO2_UV'
        hdf_path='/home/kristof/work/NDACC/HDF4_data_submission/HDF_files_UV/';
        search_str='*.no2_*.hdf';
    case 'NO2_VIS'
        hdf_path='/home/kristof/work/NDACC/HDF4_data_submission/HDF_files/';
        search_str='*.no2_*.hdf';
    case 'O3_VIS'
        hdf_path='/home/kristof/work/NDACC/HDF4_data_submission/HDF_files/';        
        search_str='*.o3_*.hdf';
    case 'O3'
        hdf_path='/home/kristof/work/bruker/';
        search_str='*.hdf';
end

tmp=dir([hdf_path search_str]);
filelist={tmp.name};

Nfiles = length(filelist);

%
for i = 1:Nfiles
    currentfile = [hdf_path filelist{i}];
    %
    % get time 
    MJDtyme = hdfread(currentfile,'DATETIME');
    % set up start and end of data index
    Nmsmts = length(MJDtyme);
    N = N + Nmsmts;
    datastart = N - Nmsmts + 1; 
    dataend = N;
    %% import data
    % time unit (DATETIME) is modified julian days 
    % (i.e. days since days since Jan. 1 2000 at 0:00:00 UT)
    tyme = datenum(2000,01,1 + MJDtyme); % in UT
    %  
    sza = double(hdfread(currentfile,'ANGLE.SOLAR_ZENITH.ASTRONOMICAL'));
    saa = double(hdfread(currentfile,'ANGLE.SOLAR_AZIMUTH'));
    z_msmt = double(hdfread(currentfile,'ALTITUDE'));

    if define_z == 1
        z = z_msmt(1,:);
        if z(1) < z(length(z))
            z = fliplr(z);
        end
    end
    %Nz = length(z);
    %P(datastart:dataend) = hdfread(currentfile,'SURFACE.PRESSURE_INDEPENDENT');
    %T(datastart:dataend) = hdfread(currentfile,'SURFACE.TEMPERATURE_INDEPENDENT');
    %%
    site.longitude = hdfread(currentfile,'LONGITUDE.INSTRUMENT');
    site.latitude = hdfread(currentfile,'LATITUDE.INSTRUMENT');
    site.altitude = hdfread(currentfile,'ALTITUDE.INSTRUMENT');
    site.name = sitename;
    %% calculate altitudes of retrieval grid points
    timevec = datevec(tyme);
    time.year(datastart:dataend) = timevec(:,1); 
    time.month(datastart:dataend) = timevec(:,2); 
    time.day(datastart:dataend) = timevec(:,3);
    time.hour(datastart:dataend) = timevec(:,4); 
    time.min(datastart:dataend) = timevec(:,5); 
    time.sec(datastart:dataend) = timevec(:,6);
    time.UTC(datastart:dataend) = zeros(length(tyme),1); % NDACC IS IN UTC
    %
    location.latitude = site.latitude;
    location.longitude = site.longitude;
    location.altitude = site.altitude;
    % FINAL MOVEMENT OF DATA FOR TOTAL RETURN 
    SZA(datastart:dataend) = sza;
    SAA(datastart:dataend) = saa;
end
fclose('all');
%
end

