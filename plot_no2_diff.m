% % % % function [ output_args ] = plot_no2_diff( input_args )
% % % %PLOT_NO2_DIFF Summary of this function goes here
% % % %   Detailed explanation goes here
% % % 
% % % %% load datasets 
% % % 
% % % % UT-GBS NO2
% % % load('/home/kristof/work/GBS/VCD_results/UT-GBS_NO2_VCD_all.mat')
% % % ut_no2=reanalysis;
% % % 
% % % % PEARL-GBS NO2 -- TEMPORARY!!
% % % load('/home/kristof/work/GBS/VCD_results/VCD_with_bad_rcd/PEARL-GBS_NO2_VCD_all.mat')
% % % reanalysis.mjd2k=mjd2k_ft(reanalysis.fd-1,reanalysis.year);
% % % reanalysis.fractional_time=reanalysis.fd-1;
% % % p_no2=reanalysis;
% % % 
% % % clearvars reanalysis
% % % 
% % % % merged GBS UV
% % % 
% % % load('/home/kristof/work/satellite_validation/GBS_NO2_UV.mat')
% % % 
% % % 
% % % %% get model NO2
% % % tmp=scale_no2_column(ut_no2.mjd2k);
% % % ut_no2.model_no2=tmp(:,1);
% % % ut_no2.scale_factor=tmp(:,2);
% % % 
% % % tmp=scale_no2_column(p_no2.mjd2k);
% % % p_no2.model_no2=tmp(:,1);
% % % p_no2.scale_factor=tmp(:,2);
% % % 
% % % tmp=scale_no2_column(gbs_no2uv.mjd2k);
% % % gbs_no2uv.model_no2=tmp(:,1);
% % % gbs_no2uv.scale_factor=tmp(:,2);

%% plot pm/am ratios

for yyyy=2003:2016
    
    figure(yyyy)
    
    subplot(311)
    hold on
    for i=unique(ut_no2.day(ut_no2.year==yyyy))'
        
        am=find(ut_no2.year==yyyy & ut_no2.day==i & ut_no2.ampm==0);
        pm=find(ut_no2.year==yyyy & ut_no2.day==i & ut_no2.ampm==1);
        
        if isempty(am+pm), continue, end
        
        diff_meas=ut_no2.mean_vcd(pm)/ut_no2.mean_vcd(am);
        diff_mod=ut_no2.model_no2(pm)/ut_no2.model_no2(am);
        
        plot(i, diff_meas, 'ro')
        plot(i, diff_mod, 'kx')
        
    end
    xlim([50,300])
    
    subplot(312)    
    hold on
    for i=unique(p_no2.day(p_no2.year==yyyy))'
        
        am=find(p_no2.year==yyyy & p_no2.day==i & p_no2.ampm==0);
        pm=find(p_no2.year==yyyy & p_no2.day==i & p_no2.ampm==1);
        
        if isempty(am+pm), continue, end
        
        diff_meas=p_no2.mean_vcd(pm)/p_no2.mean_vcd(am);
        diff_mod=p_no2.model_no2(pm)/p_no2.model_no2(am);
        
        plot(i, diff_meas, 'ro')
        plot(i, diff_mod, 'kx')
        
    end
    xlim([50,300])
    
    subplot(313)
    hold on
    for i=unique(gbs_no2uv.day(gbs_no2uv.year==yyyy))'
        
        am=find(gbs_no2uv.year==yyyy & gbs_no2uv.day==i & gbs_no2uv.ampm==0);
        pm=find(gbs_no2uv.year==yyyy & gbs_no2uv.day==i & gbs_no2uv.ampm==1);
        
        if isempty(am+pm), continue, end
        
        diff_meas=gbs_no2uv.mean_vcd(pm)/gbs_no2uv.mean_vcd(am);
        diff_mod=gbs_no2uv.model_no2(pm)/gbs_no2uv.model_no2(am);
        
        plot(i, diff_meas, 'ro')
        plot(i, diff_mod, 'kx')
        
    end
    xlim([50,300])
    

end


% end

