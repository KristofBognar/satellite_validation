function [gbs_o3,saoz_o3,bruker_o3,bruker_o3_doas,brewer_o3_ds,paris_o3,...
          ace_fts_o3,ace_mae_o3,osiris_o3,...
          ace_fts_o3_bk,ace_mae_o3_bk,osiris_o3_bk,...
          ace_fts_o3_bw,ace_mae_o3_bw,osiris_o3_bw,...
          gbs_no2,gbs_no2uv,saoz_no2,bruker_no2,...
          saoz_no2_fixRCD,saoz_no2_dailyRCD,saoz_o3_fixRCD,saoz_o3_dailyRCD,...
          saoz_no2_fixRCD2,saoz_no2_dailyRCD2,saoz_o3_fixRCD2,saoz_o3_dailyRCD2,...
          ace_fts_no2,osiris_no2,ace_fts_no2uv,osiris_no2uv,...
          ace_fts_no2_bk,osiris_no2_bk,...
          ut_o3_tmp,ut_o3_cf,saoz_o3_all,saoz_o3_cf,...
          label_gbs,label_gbs_uv,label_saoz,label_bruker,label_brewer,label_paris,label_fts,label_mae,label_osiris]...
         = load_comparison_data(alt_dmp,no2_scale,use_smooth,range_dir)
%LOAD_COMPARISON_DATA Loads data for comparison plots
%   Units are converted to be consistent
%   Data is filtered for NaNs and missing values
%   NO2 profiles are scaled to local noon
%
% kinda got out of hand...

% conversion to DU
du=2.687e16;

% select 500km radius by default
if nargin==3 
    range_dir=''; 
end

% use smoothed columns?
% use_smooth=true;

%% labels

% % label_gbs='GBS';
% % label_gbs_uv='GBS-UV';
% % label_saoz='SAOZ';
% % label_bruker='Bruker';
% % label_brewer='Brewer';
% % 
% % label_fts='ACE-FTS';
% % label_mae='MAESTRO';
% % label_osiris='OSIRIS';

label_gbs='GV';
label_gbs_uv='GU';
label_saoz='SA';
label_bruker='BK';
label_brewer='BW';
label_paris='PA';

label_fts='AF';
label_mae='AM';
label_osiris='OS';

%% load datasets 

%% CF data from Xiaoyi
load('/home/kristof/work/satellite_validation/archive_CF_project_from_Xiaoyi.mat')
ut_o3_tmp=GBS;
ut_o3_cf=GBS_CF_MERRA2_EWS;
saoz_o3_all=SAOZ;
saoz_o3_cf=SAOZ_CF_MERRA2_EWS;

ut_o3_cf.mjd2k=ft_to_mjd2k(ut_o3_cf.fd-1,ut_o3_cf.year);
ut_o3_cf.fractional_time=ut_o3_cf.fd-1;
ut_o3_cf.mean_vcd=ut_o3_cf.mean_vcd./du;
ut_o3_cf.sigma_mean_vcd=ut_o3_cf.sigma_mean_vcd./du;
ut_o3_cf.std_vcd=ut_o3_cf.std_vcd./du;

ut_o3_tmp.mjd2k=ft_to_mjd2k(ut_o3_tmp.fd-1,ut_o3_tmp.year);
ut_o3_tmp.fractional_time=ut_o3_tmp.fd-1;
ut_o3_tmp.mean_vcd=ut_o3_tmp.mean_vcd./du;
ut_o3_tmp.sigma_mean_vcd=ut_o3_tmp.sigma_mean_vcd./du;
ut_o3_tmp.std_vcd=ut_o3_tmp.std_vcd./du;

saoz_o3_cf.mjd2k=ft_to_mjd2k(saoz_o3_cf.fd-1,saoz_o3_cf.year);
saoz_o3_cf.fractional_time=saoz_o3_cf.fd-1;
saoz_o3_cf.mean_vcd=saoz_o3_cf.mean_vcd./du;
saoz_o3_cf.sigma_mean_vcd=saoz_o3_cf.sigma_mean_vcd./du;
saoz_o3_cf.std_vcd=saoz_o3_cf.std_vcd./du;

saoz_o3_all.mjd2k=ft_to_mjd2k(saoz_o3_all.fd-1,saoz_o3_all.year);
saoz_o3_all.fractional_time=saoz_o3_all.fd-1;
saoz_o3_all.mean_vcd=saoz_o3_all.mean_vcd./du;
saoz_o3_all.sigma_mean_vcd=saoz_o3_all.sigma_mean_vcd./du;
saoz_o3_all.std_vcd=saoz_o3_all.std_vcd./du;

%% merged GBS

load('/home/kristof/work/satellite_validation/GBS_NO2_UV.mat')
load('/home/kristof/work/satellite_validation/GBS_NO2.mat')
% load('/home/kristof/work/satellite_validation/GBS_O3_addRCDerr.mat')
load('/home/kristof/work/satellite_validation/GBS_O3.mat')

gbs_o3.mean_vcd=gbs_o3.mean_vcd./du;
gbs_o3.sigma_mean_vcd=gbs_o3.sigma_mean_vcd./du;
gbs_o3.std_vcd=gbs_o3.std_vcd./du;

%% SAOZ data
load('/home/kristof/work/SAOZ/saoz_o3.mat')
saoz_o3=saoz;

load('/home/kristof/work/SAOZ/saoz_no2.mat')
saoz_no2=saoz;

clearvars saoz

% SAOZ NO2 retrieved with GBS software
load('/home/kristof/work/SAOZ/VCD_results_max5_SZA/SAOZ_NO2_fixRCD_all.mat');
saoz_no2_fixRCD=data;

load('/home/kristof/work/SAOZ/VCD_results_max5_SZA/SAOZ_NO2_dailyRCD_all.mat');
saoz_no2_dailyRCD=data;

load('/home/kristof/work/SAOZ/VCD_results_max5_SZA/SAOZ_O3_fixRCD_all.mat');
data.mean_vcd=data.mean_vcd./du;
data.sigma_mean_vcd=data.sigma_mean_vcd./du;
data.std_vcd=data.std_vcd./du;
saoz_o3_fixRCD=data;

load('/home/kristof/work/SAOZ/VCD_results_max5_SZA/SAOZ_O3_dailyRCD_all.mat');
data.mean_vcd=data.mean_vcd./du;
data.sigma_mean_vcd=data.sigma_mean_vcd./du;
data.std_vcd=data.std_vcd./du;
saoz_o3_dailyRCD=data;

% SAOZ NO2 retrieved with GBS software
load('/home/kristof/work/SAOZ/VCD_results_fix_SZA/SAOZ_NO2_fixRCD_all.mat');
saoz_no2_fixRCD2=data;

load('/home/kristof/work/SAOZ/VCD_results_fix_SZA/SAOZ_NO2_dailyRCD_all.mat');
saoz_no2_dailyRCD2=data;

load('/home/kristof/work/SAOZ/VCD_results_fix_SZA/SAOZ_O3_fixRCD_all.mat');
data.mean_vcd=data.mean_vcd./du;
data.sigma_mean_vcd=data.sigma_mean_vcd./du;
data.std_vcd=data.std_vcd./du;
saoz_o3_fixRCD2=data;

load('/home/kristof/work/SAOZ/VCD_results_fix_SZA/SAOZ_O3_dailyRCD_all.mat');
data.mean_vcd=data.mean_vcd./du;
data.sigma_mean_vcd=data.sigma_mean_vcd./du;
data.std_vcd=data.std_vcd./du;
saoz_o3_dailyRCD2=data;


clearvars data

%% Bruker data
load('/home/kristof/work/bruker/bruker_o3.mat')
% load('/home/kristof/work/bruker/O3_Bavo_COMPLETE/bruker_o3.mat')
bruker_o3=bruker;
bruker_o3.tot_col=bruker_o3.tot_col./du;
bruker_o3.tot_col_err_rand=bruker_o3.tot_col_err_rand./du;
bruker_o3.tot_col_err_sys=bruker_o3.tot_col_err_sys./du;
bruker_o3.tot_col_smooth=bruker_o3.tot_col_smooth./du;
bruker_o3.part_col=bruker_o3.part_col./du;

clearvars bruker

load('/home/kristof/work/bruker/bruker_no2.mat')
bruker_no2=bruker;

%% PARIS data
load('/home/kristof/work/PARIS/PARIS_o3.mat')
paris_o3=paris;
paris_o3.tot_col=paris_o3.tot_col./du;
paris_o3.tot_col_err_rand=paris_o3.tot_col_err_rand./du;
paris_o3.tot_col_err_sys=paris_o3.tot_col_err_sys./du;
paris_o3.tot_col_smooth=paris_o3.tot_col_smooth./du;
paris_o3.part_col=paris_o3.part_col./du;

clearvars paris

%% Brewer direct sun data
load('/home/kristof/work/brewer/brewer_tables.mat')
brewer_o3_ds=brewer_ds;
brewer_o3_zs=brewer_zs;
clearvars brewer_ds brewer_zs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACE-FTS O3
load(['/home/kristof/work/satellite_validation/ACE-FTS_data/' range_dir 'ACE_v3p6_O3_table_14-52km.mat'])
ace_fts_o3=ace_fts;

%% ACE-FTS NO2
% % load('/home/kristof/work/satellite_validation/ACE-FTS_data/ACE_v3p6_NO2_table_12-32km.mat')

load(['/home/kristof/work/satellite_validation/ACE-FTS_data/' range_dir 'ACE_v3p6_NO2_table_12-40km.mat'])
ace_fts_no2=ace_fts;

clearvars ace_fts

%% ACE-MAESTRO O3
load(['/home/kristof/work/satellite_validation/ACE-MAESTRO_data/' range_dir 'ACE_MAESTRO_O3_table_14-52km.mat'])
ace_mae_o3=ace_maestro;

% load('/home/kristof/work/satellite_validation/ACE-MAESTRO_data/ACE_MAESTRO_NO2_v3p12p1_table.mat')
% ace_mae_no2=ace_maestro;

clearvars ace_maestro

%% ODIN-OSIRIS
load(['/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/' range_dir 'OSIRIS_v5p10_O3_table_14-52km.mat'])
osiris_o3=osiris;

load(['/home/kristof/work/satellite_validation/ODIN-OSIRIS_data/' range_dir 'OSIRIS_v6p00_NO2_table_12-32km.mat'])
osiris_no2=osiris;

% scale OSIRIS no2 to 40km using NDACC lut 
osiris_no2.tot_col=osiris_no2.tot_col.*osiris_no2.scaleto_40;
osiris_no2.tot_col_err=osiris_no2.tot_col_err.*osiris_no2.scaleto_40;

osiris_no2.tot_col_smooth=osiris_no2.tot_col_smooth.*osiris_no2.scaleto_40;
osiris_no2.tot_col_smooth_uv=osiris_no2.tot_col_smooth_uv.*osiris_no2.scaleto_40;
osiris_no2.tot_col_smooth_bk=osiris_no2.tot_col_smooth_bk.*osiris_no2.scaleto_40;

% osiris_no2.tot_col=osiris_no2.tot_col*1.2052; % to 40km from std atm
% osiris_no2.tot_col=osiris_no2.tot_col*1.2166; % to 42 km from std atm


clearvars osiris

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add satellite sza/saa (as seen from tangent point/height)
% ACE-FTS
date_tmp=mjd2k_to_date(ace_fts_o3.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,ace_fts_o3.lat,ace_fts_o3.lon,25*ones(size(ace_fts_o3.lon)));
sza_tmp=90-sza_tmp;

ace_fts_o3.saa=az_tmp;
ace_fts_o3.sza=sza_tmp;

date_tmp=mjd2k_to_date(ace_fts_no2.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,ace_fts_no2.lat,ace_fts_no2.lon,25*ones(size(ace_fts_no2.lon)));
sza_tmp=90-sza_tmp;

ace_fts_no2.saa=az_tmp;
ace_fts_no2.sza=sza_tmp;

% ACE-MAESTRO
date_tmp=mjd2k_to_date(ace_mae_o3.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,ace_mae_o3.lat,ace_mae_o3.lon,25*ones(size(ace_mae_o3.lon)));
sza_tmp=90-sza_tmp;

ace_mae_o3.saa=az_tmp;
ace_mae_o3.sza=sza_tmp;

% OSIRIS
date_tmp=mjd2k_to_date(osiris_o3.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,osiris_o3.lat,osiris_o3.lon,25*ones(size(osiris_o3.lon)));
sza_tmp=90-sza_tmp;

osiris_o3.saa=az_tmp;
osiris_o3.sza=sza_tmp;

date_tmp=mjd2k_to_date(osiris_no2.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,osiris_no2.lat,osiris_no2.lon,35*ones(size(osiris_no2.lon)));
sza_tmp=90-sza_tmp;

osiris_no2.saa=az_tmp;
osiris_no2.sza=sza_tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% match to DMP

[spv,temperature,theta,lat,lon]=match_DMP_DOAS(gbs_o3.fractional_time,gbs_o3.year,'O3_VIS',alt_dmp);
gbs_o3.T=temperature;
gbs_o3.spv=spv;
gbs_o3.lat_dmp=lat;
gbs_o3.lon_dmp=lon;

[spv,temperature,theta,lat,lon]=match_DMP_DOAS(gbs_no2.fractional_time,gbs_no2.year,'NO2_VIS',alt_dmp);
gbs_no2.T=temperature;
gbs_no2.spv=spv;
gbs_no2.lat_dmp=lat;
gbs_no2.lon_dmp=lon;

[spv,temperature,theta,lat,lon]=match_DMP_DOAS(gbs_no2uv.fractional_time,gbs_no2uv.year,'NO2_UV',alt_dmp);
gbs_no2uv.T=temperature;
gbs_no2uv.spv=spv;
gbs_no2uv.lat_dmp=lat;
gbs_no2uv.lon_dmp=lon;

%%%

[spv,temperature,theta,lat,lon]=match_DMP_DOAS(saoz_o3.fractional_time,saoz_o3.year,'O3_VIS',alt_dmp);
saoz_o3.T=temperature;
saoz_o3.spv=spv;
saoz_o3.lat_dmp=lat;
saoz_o3.lon_dmp=lon;

[spv,temperature,theta,lat,lon]=match_DMP_DOAS(saoz_no2.fractional_time,saoz_no2.year,'NO2_VIS',alt_dmp);
saoz_no2.T=temperature;
saoz_no2.spv=spv;
saoz_no2.lat_dmp=lat;
saoz_no2.lon_dmp=lon;

%%%

[spv,temperature,theta,lat,lon]=match_DMP_bruker(bruker_o3.fractional_time,bruker_o3.year,alt_dmp);
bruker_o3.T=temperature;
bruker_o3.spv=spv;
bruker_o3.lat_dmp=lat;
bruker_o3.lon_dmp=lon;

[spv,temperature,theta,lat,lon]=match_DMP_bruker(bruker_no2.fractional_time,bruker_no2.year,alt_dmp);
bruker_no2.T=temperature;
bruker_no2.spv=spv;
bruker_no2.lat_dmp=lat;
bruker_no2.lon_dmp=lon;

%%%

[spv,temperature,theta,lat,lon]=match_DMP_PARIS(paris_o3.fractional_time,paris_o3.year,alt_dmp);
paris_o3.T=temperature;
paris_o3.spv=spv;
paris_o3.lat_dmp=lat;
paris_o3.lon_dmp=lon;

%%%

[spv,temperature,theta,lat,lon]=match_DMP_OSIRIS(osiris_o3.fractional_time,osiris_o3.year,'O3',alt_dmp);
osiris_o3.T=temperature;
osiris_o3.spv=spv;
osiris_o3.lat_dmp=lat; %OSIRIS lat is fixed; DMPs are profiles at a single geolocation

[spv,temperature,theta,lat,lon]=match_DMP_OSIRIS(osiris_no2.fractional_time,osiris_no2.year,'NO2',alt_dmp);
osiris_no2.T=temperature;
osiris_no2.spv=spv;
osiris_no2.lat_dmp=lat;

%%%

[spv,temperature,theta,lat,lon]=match_DMP_ACE(ace_fts_o3.mjd2k,alt_dmp);
ace_fts_o3.T=temperature;
ace_fts_o3.spv=spv;
ace_fts_o3.lat_dmp=lat;

[spv,temperature,theta,lat,lon]=match_DMP_ACE(ace_fts_no2.mjd2k,alt_dmp);
ace_fts_no2.T=temperature;
ace_fts_no2.spv=spv;
ace_fts_no2.lat_dmp=lat;

[spv,temperature,theta,lat,lon]=match_DMP_ACE(ace_mae_o3.mjd2k,alt_dmp);
ace_mae_o3.T=temperature;
ace_mae_o3.spv=spv;
ace_mae_o3.lat_dmp=lat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add 30 km SZA for ground-based NO2

alt_ind=find(alt_dmp==30);

date_tmp=mjd2k_to_date(gbs_no2.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,gbs_no2.lat_dmp(:,alt_ind),...
                              gbs_no2.lon_dmp(:,alt_ind),30*ones(size(date_tmp)));
sza_tmp=90-sza_tmp;

gbs_no2.saa_30=az_tmp;
gbs_no2.sza_30=sza_tmp;

%
date_tmp=mjd2k_to_date(gbs_no2uv.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,gbs_no2uv.lat_dmp(:,alt_ind),...
                              gbs_no2uv.lon_dmp(:,alt_ind),30*ones(size(date_tmp)));
sza_tmp=90-sza_tmp;

gbs_no2uv.saa_30=az_tmp;
gbs_no2uv.sza_30=sza_tmp;

%
date_tmp=mjd2k_to_date(saoz_no2.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,saoz_no2.lat_dmp(:,alt_ind),...
                              saoz_no2.lon_dmp(:,alt_ind),30*ones(size(date_tmp)));
sza_tmp=90-sza_tmp;

saoz_no2.saa_30=az_tmp;
saoz_no2.sza_30=sza_tmp;

%
date_tmp=mjd2k_to_date(bruker_no2.mjd2k);
[az_tmp, sza_tmp] = SolarAzEl(date_tmp,bruker_no2.lat_dmp(:,alt_ind),...
                              bruker_no2.lon_dmp(:,alt_ind),30*ones(size(date_tmp)));
sza_tmp=90-sza_tmp;

bruker_no2.saa_30=az_tmp;
bruker_no2.sza_30=sza_tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% switch to smooth columns if required

bruker_o3_doas=bruker_o3;
    
ace_fts_no2uv=ace_fts_no2;
osiris_no2uv=osiris_no2;

ace_fts_o3_bk=ace_fts_o3;
ace_mae_o3_bk=ace_mae_o3;
osiris_o3_bk=osiris_o3;    

ace_fts_no2_bk=ace_fts_no2;
osiris_no2_bk=osiris_no2;    

ace_fts_o3_bw=ace_fts_o3;
ace_fts_o3_bw.tot_col=ace_fts_o3_bw.tot_col+ace_fts_o3_bw.col_below_rl;
ace_mae_o3_bw=ace_mae_o3;
ace_mae_o3_bw.tot_col=ace_mae_o3_bw.tot_col+ace_mae_o3_bw.col_below_rl;
osiris_o3_bw=osiris_o3;    
osiris_o3_bw.tot_col=osiris_o3_bw.tot_col+osiris_o3_bw.col_below_rl;

if use_smooth

    bruker_o3_doas(bruker_o3_doas.tot_col_smooth<0,:)=[];
    bruker_o3_doas.tot_col=bruker_o3_doas.tot_col_smooth./du;

    ace_fts_o3.tot_col=ace_fts_o3.tot_col_smooth;
    ace_fts_no2.tot_col=ace_fts_no2.tot_col_smooth;
    ace_mae_o3.tot_col=ace_mae_o3.tot_col_smooth;
    osiris_o3.tot_col=osiris_o3.tot_col_smooth;    
    osiris_no2.tot_col=osiris_no2.tot_col_smooth;    
    
    ace_fts_no2uv.tot_col=ace_fts_no2uv.tot_col_smooth_uv;
    osiris_no2uv.tot_col=osiris_no2uv.tot_col_smooth_uv;    

    ace_fts_no2_bk.tot_col=ace_fts_no2_bk.tot_col_smooth_bk;
    osiris_no2_bk.tot_col=osiris_no2_bk.tot_col_smooth_bk;    
    
    ace_fts_o3_bk.tot_col=ace_fts_o3.tot_col_smooth_bk;
    ace_fts_o3_bk.part_col=ace_fts_o3.part_col_smooth_bk;
    
    ace_mae_o3_bk.tot_col=ace_mae_o3.tot_col_smooth_bk;
    ace_mae_o3_bk.part_col=ace_mae_o3.part_col_smooth_bk;
    
    osiris_o3_bk.tot_col=osiris_o3.tot_col_smooth_bk;    
    osiris_o3_bk.part_col=osiris_o3.part_col_smooth_bk;    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data filters (do it individually since smoothing sometimes changes column validity)

% remove satellite NaNs (where profiles didn't cover altitude range)
% filter by max/min reasonable values if needed

o3_filter={'ace_fts_o3','ace_fts_o3_bk','ace_fts_o3_bw',...
           'ace_mae_o3','ace_mae_o3_bk','ace_mae_o3_bw',...
           'osiris_o3','osiris_o3_bk','osiris_o3_bw'};

for i=1:length(o3_filter)
    eval(['ind=find(isnan(' o3_filter{i} '.tot_col) | '...
          'isnan(' o3_filter{i} '.part_col) | '...
          o3_filter{i} '.part_col<100 | '...
          o3_filter{i} '.tot_col<100 | '...
          o3_filter{i} '.tot_col>700);'])
    eval([o3_filter{i} '(ind,:)=[];']);
end

no2_filter={'ace_fts_no2','ace_fts_no2uv','ace_fts_no2_bk',...
           'osiris_no2','osiris_no2uv','osiris_no2_bk'};

for i=1:length(no2_filter)
    eval(['ind=find(isnan(' no2_filter{i} '.tot_col) | ' no2_filter{i} '.tot_col==-9999 | '...
          no2_filter{i} '.tot_col>6e15 | ' no2_filter{i} '.tot_col<0);'])
    eval([no2_filter{i} '(ind,:)=[];']);
end

% remove NaNs from CF files
ut_o3_cf(isnan(ut_o3_cf.mean_vcd),:)=[];
ut_o3_tmp(isnan(ut_o3_tmp.mean_vcd),:)=[];
saoz_o3_cf(isnan(saoz_o3_cf.mean_vcd),:)=[];
saoz_o3_all(isnan(saoz_o3_all.mean_vcd),:)=[];

% remove CF days when std err is nan to match GBS files
ut_o3_cf(isnan(ut_o3_cf.sigma_mean_vcd),:)=[];
ut_o3_tmp(isnan(ut_o3_tmp.sigma_mean_vcd),:)=[];
saoz_o3_cf(isnan(saoz_o3_cf.sigma_mean_vcd),:)=[];
saoz_o3_all(isnan(saoz_o3_all.sigma_mean_vcd),:)=[];


% % additional filter for MAESTRO based on retrieval success percentage
% ace_mae_o3(ace_mae_o3.failed_retr_frac>0.2,:)=[];
% ace_mae_o3_bk(ace_mae_o3_bk.failed_retr_frac>0.2,:)=[];
% ace_mae_o3_bw(ace_mae_o3_bw.failed_retr_frac>0.2,:)=[];
% %ace_mae_no2(ace_mae_no2.failed_retr_frac>0.2,:)=[];

% % additional filter for MAESTRO based on reprted errors
% err=(ace_mae_o3.part_col_err./ace_mae_o3.part_col)*100;
% ace_mae_o3(err>100,:)=[];
% ace_mae_o3_bk(err>100,:)=[];
% ace_mae_o3_bw(err>100,:)=[];

% bruker NO2 filters
bruker_no2(bruker_no2.part_col<0,:)=[];


% ace_mae_no2(ace_mae_no2.tot_col<1e11,:)=[];
% ace_fts_o3(ace_fts_o3.dist>300,:)=[];
% osiris_o3(osiris_o3.dist>300,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale NO2 columns
% using interpolated albedo data:
%   interpolated sonde input gives slightly better comparison results than
%   taking nearest sonde

if no2_scale
    
%     tmp=scale_no2_column(gbs_no2.mjd2k,gbs_no2.sza_30);
    tmp=scale_no2_column(gbs_no2.mjd2k);
    gbs_no2.model_no2=tmp(:,1);
    gbs_no2.scale_factor=tmp(:,2);
    gbs_no2.mean_vcd=gbs_no2.mean_vcd.*gbs_no2.scale_factor;

%     tmp=scale_no2_column(gbs_no2uv.mjd2k,gbs_no2uv.sza_30);
    tmp=scale_no2_column(gbs_no2uv.mjd2k);
    gbs_no2uv.model_no2=tmp(:,1);
    gbs_no2uv.scale_factor=tmp(:,2);
    gbs_no2uv.mean_vcd=gbs_no2uv.mean_vcd.*gbs_no2uv.scale_factor;

%     tmp=scale_no2_column(saoz_no2.mjd2k,saoz_no2.sza_30);
    tmp=scale_no2_column(saoz_no2.mjd2k);
    saoz_no2.model_no2=tmp(:,1);
    saoz_no2.scale_factor=tmp(:,2);
    saoz_no2.mean_vcd=saoz_no2.mean_vcd.*saoz_no2.scale_factor;

    tmp=scale_no2_column(saoz_no2_fixRCD.mjd2k);
    saoz_no2_fixRCD.model_no2=tmp(:,1);
    saoz_no2_fixRCD.scale_factor=tmp(:,2);
    saoz_no2_fixRCD.mean_vcd=saoz_no2_fixRCD.mean_vcd.*saoz_no2_fixRCD.scale_factor;

    tmp=scale_no2_column(saoz_no2_dailyRCD.mjd2k);
    saoz_no2_dailyRCD.model_no2=tmp(:,1);
    saoz_no2_dailyRCD.scale_factor=tmp(:,2);
    saoz_no2_dailyRCD.mean_vcd=saoz_no2_dailyRCD.mean_vcd.*saoz_no2_dailyRCD.scale_factor;

    tmp=scale_no2_column(saoz_no2_fixRCD2.mjd2k);
    saoz_no2_fixRCD2.model_no2=tmp(:,1);
    saoz_no2_fixRCD2.scale_factor=tmp(:,2);
    saoz_no2_fixRCD2.mean_vcd=saoz_no2_fixRCD2.mean_vcd.*saoz_no2_fixRCD2.scale_factor;

    tmp=scale_no2_column(saoz_no2_dailyRCD2.mjd2k);
    saoz_no2_dailyRCD2.model_no2=tmp(:,1);
    saoz_no2_dailyRCD2.scale_factor=tmp(:,2);
    saoz_no2_dailyRCD2.mean_vcd=saoz_no2_dailyRCD2.mean_vcd.*saoz_no2_dailyRCD2.scale_factor;
    
%     tmp=scale_no2_column(bruker_no2.mjd2k,bruker_no2.sza_30);
    tmp=scale_no2_column(bruker_no2.mjd2k);
    bruker_no2.model_no2=tmp(:,1);
    bruker_no2.scale_factor=tmp(:,2);
    bruker_no2.tot_col=bruker_no2.tot_col.*bruker_no2.scale_factor;
    bruker_no2.part_col=bruker_no2.part_col.*bruker_no2.scale_factor;
    bruker_no2.part_col_doas=bruker_no2.part_col_doas.*bruker_no2.scale_factor;
    bruker_no2.tot_col_smooth=bruker_no2.tot_col_smooth.*bruker_no2.scale_factor;
    bruker_no2.tot_col_smooth_uv=bruker_no2.tot_col_smooth_uv.*bruker_no2.scale_factor;
    
    %%%
    
    tmp=scale_no2_column(ace_fts_no2.mjd2k);
    ace_fts_no2.model_no2=tmp(:,1);
    ace_fts_no2.scale_factor=tmp(:,2);
    ace_fts_no2.tot_col=ace_fts_no2.tot_col.*ace_fts_no2.scale_factor;
    ace_fts_no2.part_col=ace_fts_no2.part_col.*ace_fts_no2.scale_factor;
    ace_fts_no2.part_col_os=ace_fts_no2.part_col_os.*ace_fts_no2.scale_factor;

    tmp=scale_no2_column(ace_fts_no2uv.mjd2k);
    ace_fts_no2uv.model_no2=tmp(:,1);
    ace_fts_no2uv.scale_factor=tmp(:,2);
    ace_fts_no2uv.tot_col=ace_fts_no2uv.tot_col.*ace_fts_no2uv.scale_factor;
    ace_fts_no2uv.part_col=ace_fts_no2uv.part_col.*ace_fts_no2uv.scale_factor;

    tmp=scale_no2_column(ace_fts_no2_bk.mjd2k);
    ace_fts_no2_bk.model_no2=tmp(:,1);
    ace_fts_no2_bk.scale_factor=tmp(:,2);
    ace_fts_no2_bk.tot_col=ace_fts_no2_bk.tot_col.*ace_fts_no2_bk.scale_factor;
    ace_fts_no2_bk.part_col=ace_fts_no2_bk.part_col.*ace_fts_no2_bk.scale_factor;
    
%     tmp=scale_no2_column(ace_mae_no2.mjd2k);
%     ace_mae_no2.model_no2=tmp(:,1);
%     ace_mae_no2.scale_factor=tmp(:,2);
%     ace_mae_no2.tot_col=ace_mae_no2.tot_col.*ace_mae_no2.scale_factor;

    tmp=scale_no2_column(osiris_no2.mjd2k);
    osiris_no2.model_no2=tmp(:,1);
    osiris_no2.scale_factor=tmp(:,2);
    osiris_no2.tot_col=osiris_no2.tot_col.*osiris_no2.scale_factor;
    osiris_no2.part_col=osiris_no2.part_col.*osiris_no2.scale_factor;
    
    tmp=scale_no2_column(osiris_no2uv.mjd2k);
    osiris_no2uv.model_no2=tmp(:,1);
    osiris_no2uv.scale_factor=tmp(:,2);
    osiris_no2uv.tot_col=osiris_no2uv.tot_col.*osiris_no2uv.scale_factor;
    osiris_no2uv.part_col=osiris_no2uv.part_col.*osiris_no2uv.scale_factor;

    tmp=scale_no2_column(osiris_no2_bk.mjd2k);
    osiris_no2_bk.model_no2=tmp(:,1);
    osiris_no2_bk.scale_factor=tmp(:,2);
    osiris_no2_bk.tot_col=osiris_no2_bk.tot_col.*osiris_no2_bk.scale_factor;
    osiris_no2_bk.part_col=osiris_no2_bk.part_col.*osiris_no2_bk.scale_factor;

end



end

