% function reformat_DMPs_one_file()
%REFORMAT_DMPS read DMPs from netcdf files and save as matlab file


%% location of DMPs and saved files

% ACE
dmp_dir='/home/kristof/work/DMP/DMP_ACE_failed_GLC_filled/';
save_dir=dmp_dir;
type='ACE';


cur_dir=pwd();
cd(dmp_dir);

nc_file='ACEFTS_L2_v3p6_GLC_GEOS5MERRA2_DynEqL.nc4';

% read date info (yyyymmdd, in a cell)
% ncread doesn't reconize data type; need to read as HDF
date_cell=h5read(nc_file,'/Date');

% read time info (fractional hour)
fh=h5read(nc_file, '/Hour');

year=zeros(length(fh), 1);
mm=year;
dd=year;
HH=year;
MM=year;
SS=year;

fractional_time=year;

% Convert date and time to fractional time
for i=1:length(date_cell)
    year(i)=str2double(date_cell{i}(1:4));
    mm(i)=str2double(date_cell{i}(5:6));
    dd(i)=str2double(date_cell{i}(7:8));
    
    fh_tmp=str2double(fh{i});
    HH(i)=floor(fh_tmp);
    MM(i)=floor((fh_tmp-HH(i))*60);
    SS(i)=round((fh_tmp-HH(i)-MM(i)/60)*3600);
    
    % get fractional date (jan 1, 00:00 = 0)
    fractional_time(i)=dayofyear(year(i),mm(i),dd(i),HH(i),MM(i),SS(i))-1;

end

mjd2k=ft_to_mjd2k(fractional_time,year);

dmpstruct=struct();
dmpstruct.mjd2k=mjd2k;
dmpstruct.fractional_time=fractional_time;
dmpstruct.year=year;

% load data as arrays
dmpstruct.lat=ncread(nc_file,'Lat');
dmpstruct.lon=ncread(nc_file, 'Lon');

tmp=ncread(nc_file, 'Altitude');
dmpstruct.altitude_km=repmat(tmp,1,length(mjd2k));

dmpstruct.P=ncread(nc_file, 'Pressure');
dmpstruct.T=ncread(nc_file, 'Temperature');
dmpstruct.Theta=ncread(nc_file, 'Theta');
dmpstruct.pv=ncread(nc_file, 'PV');
dmpstruct.spv=ncread(nc_file, 'sPV');


if ~exist(save_dir,'dir'), mkdir(save_dir); end



save([save_dir type '_DMP_struct'],'dmpstruct')
    

cd(cur_dir)

% end

