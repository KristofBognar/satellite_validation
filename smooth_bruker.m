function [ bruker ] = smooth_bruker( tg, path_in, bruker, instrument )
%SMOOTH_BRUKER Summary of this function goes here
%   Detailed explanation goes here

% setup
% % tg=1;
% % load('/home/kristof/work/bruker/bruker_o3.mat')
AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';

% read bruker profiles
[alt_bk,~,~,~,prof_bk]=read_bruker_prof_avk(tg,path_in,bruker.mjd2k, instrument);


if strcmpi(tg,'o3')
    
    lowlim=.61;
    highlim=60;

    % get input for DOAS AVK LUT
    load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
    % time when we have the sonde measurement
    sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
    sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)
    % sat measurement time
    measurement_time = bruker.year + bruker.fractional_time./daysinyear(bruker.year);
    sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

    % get DOAS AVKs and aprioris
    [lut_prof,lut_avk]=read_DOAS_prof_avk(1,[bruker.year,floor(bruker.fractional_time)+1,...
                                          sonde_in],AVK_LUT_dir );

    % get avk smoothed partial column
    tot_col_smooth=integrate_smooth( alt_bk*1e5,prof_bk,...
                       lowlim*1e5,highlim*1e5,'midpoint',...
                       [0.5:59.5]*1e5, lut_prof(:,4:end), lut_avk(:,4:end), 'interp');      


    bruker.tot_col_smooth=tot_col_smooth';

elseif strcmpi(tg,'no2')
    
    lowlim=12;
    highlim=60;

    % get DOAS AVKs and aprioris (vis)
    [lut_prof,lut_avk]=read_DOAS_prof_avk(2,[bruker.year,bruker.fractional_time],...
                                          AVK_LUT_dir );

    % get avk smoothed partial column
    tot_col_smooth=integrate_smooth( alt_bk*1e5,prof_bk,...
                       lowlim*1e5,highlim*1e5,'midpoint',...
                       [0.5:59.5]*1e5, lut_prof(:,3:end), lut_avk(:,3:end), 'interp');      
    
    bruker.tot_col_smooth=tot_col_smooth';
    
    % get DOAS AVKs and aprioris (UV)
    [lut_prof,lut_avk]=read_DOAS_prof_avk(3,[bruker.year,bruker.fractional_time],...
                                          AVK_LUT_dir );

    % get avk smoothed partial column
    tot_col_smooth=integrate_smooth( alt_bk*1e5,prof_bk,...
                       lowlim*1e5,highlim*1e5,'midpoint',...
                       [0.5:59.5]*1e5, lut_prof(:,3:end), lut_avk(:,3:end), 'interp');      
    
    bruker.tot_col_smooth_uv=tot_col_smooth';
    
end

end

