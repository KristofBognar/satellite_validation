function save_bruker_prof_avk( tg, instrument )
%SAVE_BRUKER_APRIORI_AVK read aprioris and averaging kernels from bruker NDACC HDF files 
%   
% INPUT: 'O3' or 'NO2' (files must be in directories with these names)
%
% OUTPUT: saves .mat file with altitude grid (alt_km), grid spacing
% (step_km), measurement time (mjd2k), a-priori partial columns,
% retrieved partial columns and column averaging kernels (time X alt)
%

satval=0;
if nargin==1, instrument='bruker'; end

% check tg input
valid_tg={'O3','NO2','HCl','HNO3','ClONO2','HF'};
if ~any(strcmp(valid_tg,tg)), error('Select valid trace gas'); end


% location of bruker files
if satval
    filedir=['/home/kristof/work/' instrument '/' tg];
else
    
    if ~strcmpi(instrument,'bruker'), error('Set file path'), end
    
    server_path='/home/kristof/atmosp_servers/';
    save_path='/home/kristof/work/bruker/PEARL_ozone_depletion/';

    switch tg
        case 'O3'
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDF_CAMS/'];
        case 'NO2'
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDF_NO2/no2.error.fixed/'];
        otherwise
            filedir=[server_path 'data/01/eur_ndacc/HDF/HDFs_for_plotting/' lower(tg)];
    end
    
end

cur_dir=pwd();

cd(filedir)
temp = dir(['groundbased*' lower(tg) '*.hdf']); 
flist = {temp.name}; % cell array of file names

% output variables
mjd2k_meas=[];
alt_km=[];
layer_height=[];
apriori=[];
apriori_vmr=[];
avk=[];
avk_vmr=[];
dofs=[];
avk_diag=[];
sensitivity=[];
prof=[];


% load files and read data
for i=1:length(flist)

    mjd2k_tmp=hdfread(flist{i},'DATETIME');
    
    alt_tmp=hdfread(flist{i},'ALTITUDE');
    grid_bound=hdfread(flist{i},'ALTITUDE.BOUNDARIES');
    
    step_tmp=grid_bound(1,:)-grid_bound(2,:);
    
    apriori_tmp=hdfread(flist{i},[tg '.COLUMN.PARTIAL_ABSORPTION.SOLAR_APRIORI']);
    apriori_vmr_tmp=hdfread(flist{i},[tg '.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR_APRIORI']);
    avk_tmp=hdfread(flist{i},[tg '.COLUMN_ABSORPTION.SOLAR_AVK']);
    prof_tmp=hdfread(flist{i},[tg '.COLUMN.PARTIAL_ABSORPTION.SOLAR']);
    avk_full=hdfread(flist{i},[tg '.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR_AVK']);
    
    % assign variables
    mjd2k_meas=[mjd2k_meas;mjd2k_tmp'];
    alt_km=alt_tmp; % doesn't change
    layer_height=step_tmp; % doesn't change
    
    apriori=[apriori;apriori_tmp];
    apriori_vmr=[apriori_vmr;apriori_vmr_tmp*1e-6]; % convert to fraction from ppmv
    avk=[avk;avk_tmp];
    prof=[prof;prof_tmp];
    
    
    for j=1:length(mjd2k_tmp)

        % calculate degrees of freedom (full altitude range)
        dofs=[dofs; trace(squeeze(avk_full(j,:,:)))];
        
        % save  diagonal of avk
        avk_diag=[avk_diag; diag(squeeze(avk_full(j,:,:)))'];
        
        % save sensitivity (sum of avk rows at each altitude)
        sensitivity=[sensitivity; sum(squeeze(avk_full(j,:,:)),1)];
        
        % save VMR avk
        avk_vmr=cat(3,avk_vmr,fliplr(flipud(squeeze(avk_full(j,:,:)))));
    end
    
end

% flip arrays so alt increases to the right; convert to double
alt_km=double(fliplr(alt_km));
layer_height=double(fliplr(layer_height));
apriori=double(fliplr(apriori)); 
apriori_vmr=double(fliplr(apriori_vmr)); 
avk=double(fliplr(avk));
avk_vmr=double(avk_vmr);
sensitivity=double(fliplr(sensitivity));
avk_diag=double(fliplr(avk_diag));
prof=double(fliplr(prof));
dofs=double(dofs);

% convert apriori to molec/cm3 (from molec/cm2)
for i=1:length(layer_height)
    
    apriori(:,i)=apriori(:,i)/(layer_height(i)*1e5);
    prof(:,i)=prof(:,i)/(layer_height(i)*1e5);
    
end

% save results
if satval
    save(['../' instrument '_' lower(tg) '_prof_avk.mat'], 'mjd2k_meas' ,'alt_km',...
        'layer_height','apriori','apriori_vmr','avk','avk_vmr','prof','dofs',...
        'avk_diag', 'sensitivity');
else
    save([save_path instrument '_' lower(tg) '_prof_avk.mat'], 'mjd2k_meas' ,'alt_km',...
        'layer_height','apriori','apriori_vmr','avk','avk_vmr','prof','dofs',...
        'avk_diag', 'sensitivity');
    
end

cd(cur_dir)

end

