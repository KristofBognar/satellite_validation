function [ alt_out, layer_height_out, ap_out, avk_out, prof_out, dof_out, avk_diag_out ] =...
         read_bruker_prof_avk( tg, path_in, mjd2k_in, instrument )
%READ_BRUKER_APRIORI_AVK read aprioris and averaging kernels from bruker NDACC HDF files 
%   
% INPUT:
%   tg: 1 for ozone, 2 for NO2
%   mjd2k_in: time
%
% OUTPUT: profiles corresponding to nearest bruker measurement
%   alt_out: altitude grid in km (layer centers)
%   layer_height_out: layer height in km (uneven)
%   ap_out: bruker apriori profile in molec/cm3
%   avk_out: column averaging kernel
%   prof_out: retrieved profile in molec/cm3

if nargin==2, instrument='bruker'; end

% load avk file
load([path_in instrument '_' lower(tg) '_prof_avk.mat'])


alt_out=alt_km;
layer_height_out=layer_height;
ap_out=NaN(length(mjd2k_in),length(alt_out));
avk_out=NaN(length(mjd2k_in),length(alt_out));
avk_diag_out=NaN(length(mjd2k_in),length(alt_out));
prof_out=NaN(length(mjd2k_in),length(alt_out));
dof_out=NaN(length(mjd2k_in),1);

% find closest bruker measurement
for i=1:length(mjd2k_in)

    [~,ind]=sort(abs(mjd2k_meas-mjd2k_in(i)));
    ind=ind(1);
    
    ap_out(i,:)=apriori(ind,:);
    avk_out(i,:)=avk(ind,:);
    avk_diag_out(i,:)=avk_diag(ind,:);
    prof_out(i,:)=prof(ind,:);
    dof_out(i,:)=dofs(ind);

end

end

