% function reformat_GLCs()
%REFORMAT_GLCS read ACE geolocation files from netcdf files and save as matlab file


%% location of DMPs and saved files

dmp_dir='/home/kristof/work/DMP/DMP_ACE/';
save_dir=dmp_dir;

cur_dir=pwd();
cd(dmp_dir);

nc_file='ACEFTS_L2_v3p6_GLC.nc';

% read date and time info
year=double(ncread(nc_file,'year'));
mm=double(ncread(nc_file,'month'));
dd=double(ncread(nc_file,'day'));
fh=double(ncread(nc_file, 'hour'));

fractional_time=zeros(length(fh), 1);

% Convert date and time to fractional time
for i=1:length(fh)
    HH=floor(fh(i));
    MM=floor((fh(i)-HH)*60);
    SS=round((fh(i)-HH-MM/60)*3600);
    
    % get fractional date (jan 1, 00:00 = 0)
    fractional_time(i)=dayofyear(year(i),mm(i),dd(i),HH,MM,SS)-1;

end

mjd2k=ft_to_mjd2k(fractional_time,year);

glcstruct=struct();
glcstruct.mjd2k=mjd2k;
glcstruct.fractional_time=fractional_time;
glcstruct.year=year;

% load data as arrays
glcstruct.lat=ncread(nc_file,'GLC_latitude');
glcstruct.lon=ncread(nc_file, 'GLC_longitude');

glcstruct.lat_30km=ncread(nc_file,'latitude');
glcstruct.lon_30km=ncread(nc_file, 'longitude');

glcstruct.QF=ncread(nc_file, 'quality_flag');

tmp=ncread(nc_file, 'altitude');
glcstruct.altitude_km=repmat(tmp,1,length(mjd2k));

glcstruct.P=ncread(nc_file, 'pressure');
glcstruct.T=ncread(nc_file, 'temperature');


if ~exist(save_dir,'dir'), mkdir(save_dir); end

save([save_dir 'ACE_GLC_struct'],'glcstruct')
    

cd(cur_dir)

% end

