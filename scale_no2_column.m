function [ out ] = scale_no2_column( mjd2k, sza )
%SCALE_NO2_COLUMN scale NO2 columns to local noon using box model output
%

% in: mjd2k; non-local SZA optional
% out: no2, scale factor

%%% used for ozone depletion paper
load('/home/kristof/work/models/NO2_box_model/no2_cols_10-60_km_feb-may_albedo_interp.mat')

%%% used for satval paper
%%%load('/home/kristof/work/models/NO2_box_model/no2_cols_12-40_km_albedo_interp.mat')

% load('/home/kristof/work/models/NO2_box_model/no2_cols_14-40_km_albedo_interp.mat')
% load('/home/kristof/work/models/NO2_box_model/no2_cols_14-40_km_nearest_sonde_albedo_interp.mat')


out=zeros(length(mjd2k),2);

for i=1:length(mjd2k)
    
    % find corresponding time
    [min_val,ind_match]=min(abs(model_no2.mjd2k-mjd2k(i)));
    
    if nargin==1
        % only time is provided, use that to get no2 scale factor
        match=model_no2.no2(ind_match);
    else
        % SZA is also provided -- assumed to be SZA at a point other than
        % the instrument location (e.g. 30 km line of sight point for ground-based
        % instruments)
        % match SZA for same twiligt indicated by time
        ind_tmp=find(model_no2.year==model_no2.year(ind_match) & ...
              model_no2.day==model_no2.day(ind_match) & ...
              model_no2.ampm==model_no2.ampm(ind_match));

        [~,ind_match_sza]=min(abs(model_no2.sza(ind_tmp)-sza(i)));
        
        % match to index in entire array
        ind_match_sza=ind_tmp(1)-1+ind_match_sza;
        
        match=model_no2.no2(ind_match_sza);
        
    end
    
    % find corresponding noon
    ind_noon=find(model_no2.year==model_no2.year(ind_match) & ...
                  model_no2.day==model_no2.day(ind_match) & model_no2.ampm==1);

    noon=model_no2.no2(ind_noon(1));
              
    ratio=noon/match;
    
    % double check if model data actually exists for given time
    if min_val<1 % data within 1 day
        out(i,:)=[match,ratio];
    else
        out(i,:)=[NaN,NaN];
    end
end




end

