function refresh_data(range_dir)
%REFRESH_DATA updates databases for instrument data
%   

if nargin==0, range_dir=''; end

% for use_smooth=0:1
for use_smooth=0
    disp('Skipping AVK smoothing step')

    % get all data (filtered, no2 scaled, avk smoothed)
    no2_scale=true;
    alt_dmp=[13.5,14,14.5,17.5,18,18.5,19.5,20,20.5,21.5,22,22.5,24,25,26,30,36];
    
    [gbs_o3,saoz_o3,bruker_o3,bruker_o3_doas,brewer_o3_ds,paris_o3,...
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
    = load_comparison_data(alt_dmp,no2_scale,use_smooth,range_dir);

    if ~isempty(range_dir)
        range_save=['_' range_dir(1:end-1)]; 
    else
        range_save=range_dir;
    end

    
    if use_smooth
        % save data
        save(['/home/kristof/work/satellite_validation/all_data_smooth' range_save '.mat'])
    else
        % save data
        save(['/home/kristof/work/satellite_validation/all_data_nosmooth' range_save '.mat'])
        
        % save tc results
        if isempty(range_dir)
            % get triple-collocation results for not smoothed case
            [tc_o3_rt,tc_o3_rmse,tc_o3_snr,tc_o3_N,tc_o3_mean_all,tc_o3_std_all]=get_tc_stats_o3(...
                  gbs_o3,saoz_o3,bruker_o3,paris_o3,brewer_o3_ds,osiris_o3,ace_fts_o3,ace_mae_o3);

            [tc_no2_rt,tc_no2_rmse,tc_no2_snr,tc_no2_N,tc_no2_mean_all,tc_no2_std_all]=get_tc_stats_no2(...
                  gbs_no2,gbs_no2uv,saoz_no2,bruker_no2,osiris_no2,ace_fts_no2);
        
        save(['/home/kristof/work/satellite_validation/all_data_tc' range_save '.mat'],...
            'tc_o3_rt','tc_o3_rmse','tc_o3_snr','tc_o3_N','tc_o3_mean_all','tc_o3_std_all',...
            'tc_no2_rt','tc_no2_rmse','tc_no2_snr','tc_no2_N','tc_no2_mean_all','tc_no2_std_all')
              
        end
        
    end
        
end

end

