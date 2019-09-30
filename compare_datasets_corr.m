function compare_datasets_corr()
% function to compare various ground-based and satellite datasets and produce
% correlation plots with best fit lines
%
% Kristof Bognar, December 2017 

% what to compare
o3_tc=0;
o3_tc_cf=0;
o3_gb=0;
no2_tc=0;
no2_gb=0;

% needs to be run separately from gb comparisons
o3_sat=0;
% no2_sat=0;

test=0;

% day range for seasonal plots where the end of 'spring' is the last day
% when the sun dips below the horizon
% day_range=[40.25,105.25];
% day_range=[106.25,239.25];
% day_range=[240.25,365.25];
% year_range=[2000,2011];

% load data
if o3_sat 
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
else
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
end

if exist(data_file,'file'), load(data_file); end

%% filter by time
if exist('day_range','var')
    % GBS and saoz
    ind=find(gbs_o3.fractional_time<day_range(1)-1 | gbs_o3.fractional_time>day_range(2));
    gbs_o3(ind,:)=[];

    ind=find(gbs_no2.fractional_time<day_range(1)-1 | gbs_no2.fractional_time>day_range(2));
    gbs_no2(ind,:)=[];

    ind=find(gbs_no2uv.fractional_time<day_range(1)-1 | gbs_no2uv.fractional_time>day_range(2));
    gbs_no2uv(ind,:)=[];

    ind=find(saoz_o3.fractional_time<day_range(1)-1 | saoz_o3.fractional_time>day_range(2));
    saoz_o3(ind,:)=[];

    ind=find(saoz_no2.fractional_time<day_range(1)-1 | saoz_no2.fractional_time>day_range(2));
    saoz_no2(ind,:)=[];

    ind=find(saoz_no2_fixRCD.fractional_time<day_range(1)-1 |...
             saoz_no2_fixRCD.fractional_time>day_range(2));
    saoz_no2_fixRCD(ind,:)=[];
    
    % bruker and brewer
    ind=find(brewer_o3_ds.fractional_time<day_range(1)-1 | brewer_o3_ds.fractional_time>day_range(2));
    brewer_o3_ds(ind,:)=[];

    ind=find(bruker_o3.fractional_time<day_range(1)-1 | bruker_o3.fractional_time>day_range(2));
    bruker_o3(ind,:)=[];

end

if exist('year_range','var')
    % GBS and saoz
    ind=find(gbs_o3.year<year_range(1)-1 | gbs_o3.year>year_range(2)-1);
    gbs_o3(ind,:)=[];

    ind=find(gbs_no2.year<year_range(1)-1 | gbs_no2.year>year_range(2)-1);
    gbs_no2(ind,:)=[];

    ind=find(gbs_no2uv.year<year_range(1)-1 | gbs_no2uv.year>year_range(2));
    gbs_no2uv(ind,:)=[];

    ind=find(saoz_o3.year<year_range(1)-1 | saoz_o3.year>year_range(2)-1);
    saoz_o3(ind,:)=[];

    ind=find(saoz_no2.year<year_range(1)-1 | saoz_no2.year>year_range(2)-1);
    saoz_no2(ind,:)=[];

    ind=find(saoz_no2_fixRCD.year<year_range(1)-1 | saoz_no2_fixRCD.year>year_range(2)-1);
    saoz_no2_fixRCD(ind,:)=[];

    % bruker and brewer
    ind=find(brewer_o3_ds.year<year_range(1)-1 | brewer_o3_ds.year>year_range(2)-1);
    brewer_o3_ds(ind,:)=[];

    ind=find(bruker_o3.year<year_range(1)-1 | bruker_o3.year>year_range(2)-1);
    bruker_o3(ind,:)=[];

end


%% Correlation plots
%% O3 total columns 
if o3_tc
    
    disp('O3 TC')
    
    figure()
    
    % w=900 for 2x2, 1200 for 2x3
%     set(gcf, 'Position', [100, 100, 900, 650]);
    
    plotrow=3;
    plotcol=4;
    ii=1;
    fprintf('\n')

    % OSIRIS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs OSIRIS')
    tmp_time(gbs_o3,osiris_o3,label_gbs,label_osiris,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ O3 vs OSIRIS O3')
    tmp_time(saoz_o3,osiris_o3,label_saoz,label_osiris,'');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Bruker vs OSIRIS')
    tmp_time(bruker_o3,osiris_o3_bk,label_bruker,label_osiris,'');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Brewer vs OSIRIS O3')
    tmp_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'bottom');
    fprintf('\n')

    % ACE-FTS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-FTS all')
    tmp_twilight(gbs_o3,ace_fts_o3,label_gbs,label_fts,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-FTS all')
    tmp_twilight(saoz_o3,ace_fts_o3,label_saoz,label_fts,'');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Bruker vs ACE-FTS all')
    tmp_time(bruker_o3,ace_fts_o3_bk,label_bruker,label_fts,'');
    fprintf('\n')
    
    ii=ii+1;


    % ACE-MAESTRO vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-MAESTRO')
    tmp_twilight(gbs_o3,ace_mae_o3,label_gbs,label_mae,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-MAESTRO')
    tmp_twilight(saoz_o3,ace_mae_o3,label_saoz,label_mae,'bottom');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Bruker vs ACE-MAESTRO')
    tmp_time(bruker_o3,ace_mae_o3_bk,label_bruker,label_mae,'bottom');
    fprintf('\n')

% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('Brewer vs ACE-FTS O3 all')
% % %     tmp_time(brewer_o3_ds,ace_fts_o3_bw,label_brewer,label_fts,'bottomleft');
% % %     fprintf('\n')
% % %     
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('Brewer vs ACE-FTS O3 all')
% % %     tmp_time(brewer_o3_ds,ace_mae_o3_bw,label_brewer,label_mae,'bottomleft');
% % %     fprintf('\n')

end

%% O3 total columns with cloud filter for GB measurements
if o3_tc_cf
    
    disp('O3 TC CF')
    
    figure()
    plotrow=3;
    plotcol=4;
    ii=1;
    fprintf('\n')

%     ut_o3_tmp=ut_o3(ut_o3.year>=2010,:);
    
    % OSIRIS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs OSIRIS')
    tmp_time(ut_o3_tmp,osiris_o3,label_gbs,label_osiris,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ O3 vs OSIRIS O3')
    tmp_time(saoz_o3_all,osiris_o3,label_saoz,label_osiris,'');
    fprintf('\n')
    
    % OSIRIS vs GB CF
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs OSIRIS')
    tmp_time(ut_o3_cf,osiris_o3,label_gbs,label_osiris,'');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ O3 vs OSIRIS O3')
    tmp_time(saoz_o3_cf,osiris_o3,label_saoz,label_osiris,'');
    fprintf('\n')

    % ACE-FTS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-FTS all')
    tmp_twilight(ut_o3_tmp,ace_fts_o3,label_gbs,label_fts,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-FTS all')
    tmp_twilight(saoz_o3_all,ace_fts_o3,label_saoz,label_fts,'');
    fprintf('\n')

    % ACE-FTS vs GB CF
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-FTS all')
    tmp_twilight(ut_o3_cf,ace_fts_o3,label_gbs,label_fts,'');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-FTS all')
    tmp_twilight(saoz_o3_cf,ace_fts_o3,label_saoz,label_fts,'');
    fprintf('\n')

    % ACE-MAESTRO vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-MAESTRO')
    tmp_twilight(ut_o3_tmp,ace_mae_o3,label_gbs,label_mae,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-MAESTRO')
    tmp_twilight(saoz_o3_all,ace_mae_o3,label_saoz,label_mae,'bottom');
    fprintf('\n')

    % ACE-MAESTRO vs GB CF
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-MAESTRO')
    tmp_twilight(ut_o3_cf,ace_mae_o3,[label_gbs '_C_F'],label_mae,'bottom');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-MAESTRO')
    tmp_twilight(saoz_o3_cf,ace_mae_o3,[label_saoz '_C_F'],label_mae,'bottom');
    fprintf('\n')

end

%% O3 satellite partial columns
if o3_sat
    
    disp('O3 SAT')
    
    figure()
    set(gcf, 'Position', [100, 100, 1200, 650]);    
    
    plotrow=2;
    plotcol=3;
    ii=1;
    fprintf('\n')
    
    % OSIRIS vs others
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-FTS vs OSIRIS')
    tmp_time(ace_fts_o3,osiris_o3,label_fts,label_osiris,'bottomleft', 12, true);
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-MAESTRO vs OSIRIS')
    tmp_time(ace_mae_o3,osiris_o3,label_mae,label_osiris,'bottomleft', 12, true);
    fprintf('\n')
    
%     ii=ii+1;
    
    % ACE-FTS vs others
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-MAESTRO vs ACE-FTS')
    tmp_twilight(ace_mae_o3,ace_fts_o3,label_mae,label_fts,'bottomleft', true);
    fprintf('\n')

%     ii=ii+1;
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-FTS vs OSIRIS')
%     tmp_time(osiris_o3,ace_fts_o3,label_osiris,label_fts,'bottomleft', 12, true);
%     fprintf('\n')
% 
%     % ACE-MAESTRO vs others
%     ii=ii+1;
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs ACE-FTS')
%     tmp_twilight(ace_fts_o3,ace_mae_o3,label_fts,label_mae,'bottomleft', true);
%     fprintf('\n')
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs OSIRIS')
%     tmp_time(osiris_o3,ace_mae_o3,label_osiris,label_mae,'bottomleft', 12, true);
%     fprintf('\n')
    
    % include NO2 plot here, since it's only a single pair
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-FTS vs OSIRIS')
    ace_fts_no2.part_col=ace_fts_no2.part_col_os;
    tmp_time(ace_fts_no2,osiris_no2,label_fts,label_osiris,'bottomleft', 12, true);
    fprintf('\n')

    
    
end

%% O3 from ground-based instruments
if o3_gb
    
    disp('O3 GB')
    
    figure()
% %     plotrow=5;
% %     plotcol=5;
    ii=1;
    fprintf('\n')
    
    fig_ax = tight_subplot(4,4,[0.05,0.03],[0.1,0.1],[0.1,0.1]);
    axes(fig_ax(16))
    delete(gca)
    
    % brewer on x-axis
    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,saoz_o3_fixRCD,label_brewer,[label_saoz '_{fix}'],'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,saoz_o3_dailyRCD,label_brewer,[label_saoz '_{day}'],'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'bottomleft',12,false);
    ii=ii+1;

    % bruker on x-axis
    axes(fig_ax(ii))
    tmp_time(bruker_o3,gbs_o3,label_bruker,label_gbs,'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(bruker_o3,saoz_o3,label_bruker,label_saoz,'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(bruker_o3,saoz_o3_fixRCD,label_bruker,[label_saoz '_{fix}'],'bottomleft',12,false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(bruker_o3,saoz_o3_dailyRCD,label_bruker,[label_saoz '_{day}'],'bottomleft',12,false);
    ii=ii+1;
    
    % saoz on x-axis
    axes(fig_ax(ii))
    tmp_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_twilight(saoz_o3,saoz_o3_fixRCD,label_saoz,[label_saoz '_{fix}'],'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_twilight(saoz_o3,saoz_o3_dailyRCD,label_saoz,[label_saoz '_{day}'],'bottomleft',false);
    ii=ii+1;
    
    % GBS on x-axis
    axes(fig_ax(ii))
    tmp_twilight(saoz_o3_fixRCD,gbs_o3,[label_saoz '_{fix}'],label_gbs,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_twilight(saoz_o3_dailyRCD,gbs_o3,[label_saoz '_{day}'],label_gbs,'bottomleft',false);
    ii=ii+1;
    
    % saoz fix vs daily
    axes(fig_ax(ii))
    tmp_twilight(saoz_o3_fixRCD,saoz_o3_dailyRCD,[label_saoz '_{fix}'],[label_saoz '_{day}'],'bottomleft',false);
    ii=ii+1;

% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs GBS O3')
% %     tmp_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs GBS O3')
% %     tmp_time(bruker_o3_doas,gbs_o3,label_bruker,label_gbs,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 all vs GBS O3')
% %     tmp_twilight(saoz_o3_all,gbs_o3,'SA2',label_gbs,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 vs GBS O3')
% %     tmp_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottomleft');
% %     fprintf('\n')
% %     
% %     ii=ii+1;
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs SAOZ O3')
% %     tmp_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs SAOZ O3')
% %     tmp_time(bruker_o3_doas,saoz_o3,label_bruker,label_saoz,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 all vs SAOZ O3')
% %     tmp_twilight(saoz_o3_all,saoz_o3,'SA2',label_saoz,'bottomleft');
% %     fprintf('\n')
% %     
% %     ii=ii+1;
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 vs GBS O3')
% %     tmp_twilight(gbs_o3,saoz_o3,label_gbs,label_saoz,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs SAOZ O3 all')
% %     tmp_time(brewer_o3_ds,saoz_o3_all,label_brewer,'SA2','bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs SAOZ O3 all')
% %     tmp_time(bruker_o3_doas,saoz_o3_all,label_bruker,'SA2','bottomleft');
% %     fprintf('\n')
% %     
% %     ii=ii+1;
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 vs SAOZ O3 all')
% %     tmp_twilight(saoz_o3,saoz_o3_all,label_saoz,'SA2','bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ O3 all vs GBS O3')
% %     tmp_twilight(gbs_o3,saoz_o3_all,label_gbs,'SA2','bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs Brewer O3')
% %     tmp_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'bottomleft');
% %     fprintf('\n')
% %     
% %     ii=ii+1;
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs SAOZ O3 all')
% %     tmp_time(saoz_o3_all,bruker_o3_doas,'SA2',label_bruker,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs SAOZ O3')
% %     tmp_time(saoz_o3,bruker_o3_doas,label_saoz,label_bruker,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs GBS O3')
% %     tmp_time(gbs_o3,bruker_o3_doas,label_gbs,label_bruker,'bottomleft');
% %     fprintf('\n')
% % 
% %     ii=ii+1;
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker O3 vs Brewer O3')
% %     tmp_time(bruker_o3,brewer_o3_ds,label_bruker,label_brewer,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs SAOZ O3 all')
% %     tmp_time(saoz_o3_all,brewer_o3_ds,'SA2',label_brewer,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs SAOZ O3')
% %     tmp_time(saoz_o3,brewer_o3_ds,label_saoz,label_brewer,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Brewer O3 vs GBS O3')
% %     tmp_time(gbs_o3,brewer_o3_ds,label_gbs,label_brewer,'bottomleft');
% %     fprintf('\n')

end

    
%% NO2
if no2_tc
    
    disp('NO2 TC')
    
    figure()
    
    % w=900 for 2x2; 1200 for 2x3
%     set(gcf, 'Position', [100, 100, 1200, 650]);
%     set(gcf, 'Position', [100, 100, 900, 700]);
    
    plotrow=3;
    plotcol=3;
    ii=1;
    fprintf('\n')

    
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('Bruker vs OSIRIS')
%     tmp_time(bruker_no2,osiris_no2_bk,label_bruker,label_osiris,'');
%     fprintf('\n')

    % ACE-FTS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-FTS all')
    tmp_twilight(gbs_no2,ace_fts_no2,label_gbs,label_fts,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GU vs ACE-FTS all')
    tmp_twilight(gbs_no2uv,ace_fts_no2uv,label_gbs_uv,label_fts,'');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-FTS all')
    tmp_twilight(saoz_no2,ace_fts_no2,label_saoz,label_fts,'');
    fprintf('\n')

    % OSIRIS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs OSIRIS')
    tmp_time(gbs_no2,osiris_no2,label_gbs,label_osiris,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GU vs OSIRIS')
    tmp_time(gbs_no2uv,osiris_no2uv,label_gbs_uv,label_osiris,'bottom');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs OSIRIS')
    tmp_time(saoz_no2,osiris_no2,label_saoz,label_osiris,'bottom');
    fprintf('\n')

    % OSIRIS vs SAOZ all year
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ fixRCD vs OSIRIS')
    tmp_time(saoz_no2_fixRCD,osiris_no2,'SAOZ_{fix}',label_osiris,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ dailyRCD vs OSIRIS O3')
    tmp_time(saoz_no2_dailyRCD,osiris_no2,'SAOZ_{day}',label_osiris,'bottom');
    fprintf('\n')
    
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('Bruker vs ACE-FTS all')
%     tmp_time(bruker_no2,ace_fts_no2,label_bruker,label_fts,'');
%     fprintf('\n')

% % %     % ACE-MAESTRO vs GB
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('UT-GBS vs ACE-MAESTRO')
% % %     tmp_twilight(gbs_no2,ace_mae_no2,label_gbs,label_mae,'bottomleft');
% % %     fprintf('\n')
% % % 
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('GU vs ACE-MAESTRO')
% % %     tmp_twilight(gbs_no2uv,ace_mae_no2,label_gbs_uv,label_mae,'bottomleft');
% % %     fprintf('\n')
% % %     
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('SAOZ vs ACE-MAESTRO')
% % %     tmp_twilight(saoz_no2,ace_mae_no2,label_saoz,label_mae,'bottom');
% % %     fprintf('\n')

%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('Bruker vs ACE-MAESTRO')
%     tmp_time(bruker_no2,ace_mae_no2,label_bruker,label_mae,'bottom');
%     fprintf('\n')



    fprintf('\n')
end

%% NO2 from ground-based instruments
if no2_gb
    
    disp('NO2 GB')
    
    figure()
    plotrow=4;
    plotcol=4;
    ii=1;
    fprintf('\n')

    bruker_no2.tot_col=bruker_no2.part_col;    
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT NO2 vs GBS NO2 UV')
    tmp_twilight(gbs_no2,gbs_no2uv,label_gbs,label_gbs_uv,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ NO2 vs GBS NO2 UV')
    tmp_twilight(saoz_no2,gbs_no2uv,label_saoz,label_gbs_uv,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ fixRCD vs GBS NO2 UV')
    tmp_twilight(saoz_no2_fixRCD,gbs_no2uv,'SAOZ_{fix}',label_gbs_uv,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ dailyRCD vs GBS NO2 UV')
    tmp_twilight(saoz_no2_dailyRCD,gbs_no2uv,'SAOZ_{day}',label_gbs_uv,'bottomleft');
    fprintf('\n')
    
    %%%

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('BK vs GBS NO2 UV')
    tmp_time(bruker_no2,gbs_no2uv,label_bruker,label_gbs_uv,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ NO2 vs GBS NO2')
    tmp_twilight(saoz_no2,gbs_no2,label_saoz,label_gbs,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ fixRCD vs GBS NO2')
    tmp_twilight(saoz_no2_fixRCD,gbs_no2,'SAOZ_{fix}',label_gbs,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ dailyRCD vs GBS NO2')
    tmp_twilight(saoz_no2_dailyRCD,gbs_no2,'SAOZ_{day}',label_gbs,'bottomleft');
    fprintf('\n')
    
    %%%

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('BK vs GBS NO2')
    tmp_time(bruker_no2,gbs_no2,label_bruker,label_gbs,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ fixRCD vs SAOZ dailyRCD')
    tmp_twilight(saoz_no2_fixRCD,saoz_no2_dailyRCD,'SAOZ_{fix}','SAOZ_{day}','bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs SAOZ fixRCD')
    tmp_twilight(saoz_no2,saoz_no2_fixRCD,label_saoz,'SAOZ_{fix}','bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs SAOZ dailyRCD')
    tmp_twilight(saoz_no2,saoz_no2_dailyRCD,label_saoz,'SAOZ_{day}','bottomleft');
    fprintf('\n')
    
    %%%
    
    ii=ii+1;
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('BK vs SAOZ')
    tmp_time(bruker_no2,saoz_no2,label_bruker,label_saoz,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('BK vs SAOZ fixRCD')
    tmp_time(bruker_no2,saoz_no2_fixRCD,label_bruker,'SAOZ_{fix}','bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('BK vs SAOZ dailyRCD')
    tmp_time(bruker_no2,saoz_no2_dailyRCD,label_bruker,'SAOZ_{day}','bottomleft');
    fprintf('\n')
    

end

if test
    
% %     figure
% %     plotrow=3;
% %     plotcol=1;
% %     ii=1;
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('UT-GBS vs ACE-FTS')
% %     tmp_twilight(ace_fts_o3,gbs_o3,'ACE',label_gbs,'', false, true);
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ vs ACE-FTS all')
% %     tmp_twilight(ace_fts_o3,saoz_o3,'ACE',label_saoz,'', false, true);
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('Bruker vs ACE-FTS all')
% %     tmp_time(ace_fts_o3,bruker_o3,'ACE',label_bruker,'', 12, false, true);
% %     fprintf('\n')

figure

% fts_rel_err=(ace_fts_o3.tot_col_err./ace_fts_o3.tot_col)*100;
% gbs_rel_err=(sqrt(gbs_o3.sigma_mean_vcd.^2+gbs_o3.std_vcd.^2)./gbs_o3.mean_vcd)*100;


% [ x1, x2, e1, e2, ~ ] = find_coincidences_time( osiris_o3, brewer_o3_ds );
[ x1, x2, e1, e2, ~ ] = find_coincidences_twilight( ace_fts_o3, gbs_o3 );

rel_err1=(e1./x1)*100;
rel_err2=(e2./x2)*100;
tot_rel_diff=(2*(x1-x2)./(x1+x2))*100;

bins=[0.5:0.5:5, 6:10, 12:2:30];

[counts,centers]=hist(tot_rel_diff,[fliplr(bins)*-1,0,bins]);
counts=counts/sum(counts);
counts(counts<1e-3)=0;
plot(centers,counts,'k-', 'linewidth',2), hold on

[counts,centers]=hist(rel_err1,bins);
counts=counts/sum(counts);
counts(counts<1e-3)=0;
plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'b-'), hold on

[counts,centers]=hist(rel_err2,bins);
counts=counts/sum(counts);
counts(counts<1e-3)=0;
plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'r-'), hold on

set(gca, 'YScale', 'log')


end

end

%% local functions to avoid copy-pasting code
function tmp_twilight(table1, table2, instr1, instr2, loc_str, partcol)

    if nargin<6, partcol=false; end

    [ x1, x2, e1, e2, times, DMP] =...
        find_coincidences_twilight( table1, table2, partcol );

    T2=DMP.T2(:,8);
    spv1=DMP.spv1(:,8);
    lat1=DMP.lat1(:,end);
    lat2=DMP.lat2;
    
    ind=[];
%     ind=find(spv1>1.2e-4 & spv1<1.6e-4);
%     ind=find(abs(lat1-lat2)>2);
    x1(ind)=[];
    x2(ind)=[];
    e1(ind)=[];
    e2(ind)=[];

    corr_plot([x1,e1],[x2,e2],instr1, instr2, loc_str, partcol);
    
    print_diffs(x1,x2)
    
%     [~,u1,u2]=meas_uncertainty(x1,x2);
%     
%     fprintf('Uncertainty 1: %3.5g\n', 100*u1/mean(x1))
%     fprintf('Uncertainty 2: %3.5g\n', 100*u2/mean(x2))
    
    
% %     rel_err1=(e1./x1)*100;
% %     rel_err2=(e2./x2)*100;
% %     tot_rel_diff=(2*(x1-x2)./(x1+x2))*100;
% % 
% %     bins=[0.5:0.5:5, 6:10, 12:2:30];
% % 
% %     [counts,centers]=hist(tot_rel_diff,[fliplr(bins)*-1,0,bins]);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot(centers,counts,'k-', 'linewidth',2), hold on
% % 
% %     [counts,centers]=hist(rel_err1,bins);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'b-'), hold on
% % 
% %     [counts,centers]=hist(rel_err2,bins);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'r-'), hold on
% % 
% %     set(gca, 'YScale', 'log')
% %     
% %     xlabel('Rel. diff. (%)')
% %     ylabel([instr2 '-' instr1])
% %     xlim([-30,30])


    
    
end

function tmp_time(table1, table2, instr1, instr2, loc_str, dt, partcol)

    if nargin==5, partcol=false; dt=12; end
    
    [ x1, x2, e1, e2, times,~, DMP ] =...
        find_coincidences_time( table1, table2, dt, partcol );

    T2=DMP.T2(:,8);
    spv1=DMP.spv1(:,8);
    lat1=DMP.lat1(:,end);
    lat2=DMP.lat2;

    ind=[];
%     ind=find(spv1>1.2e-4 & spv1<1.6e-4);
%     ind=find(abs(lat1-lat2)>2);
    x1(ind)=[];
    x2(ind)=[];
    e1(ind)=[];
    e2(ind)=[];

    corr_plot([x1,e1],[x2,e2],instr1, instr2, loc_str, partcol);
    
    print_diffs(x1,x2);
    
%     disp(['Mean abs. diff.: ' num2str(abs_diff)])
%     disp(['Mean rel. diff.: ' num2str(rel_diff) '%'])

% %     rel_err1=(e1./x1)*100;
% %     rel_err2=(e2./x2)*100;
% %     tot_rel_diff=(2*(x1-x2)./(x1+x2))*100;
% % 
% %     bins=[0.5:0.5:5, 6:10, 12:2:30];
% % 
% %     [counts,centers]=hist(tot_rel_diff,[fliplr(bins)*-1,0,bins]);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot(centers,counts,'k-', 'linewidth',2), hold on
% % 
% %     [counts,centers]=hist(rel_err1,bins);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'b-'), hold on
% % 
% %     [counts,centers]=hist(rel_err2,bins);
% %     counts=counts/sum(counts);
% %     counts(counts<1e-3)=0;
% %     plot([fliplr(centers)*-1 centers],[fliplr(counts) counts],'r-'), hold on
% % 
% %     set(gca, 'YScale', 'log')
% %     
% %     xlabel('Rel. diff. (%)')
% %     ylabel([instr2 '-' instr1])
% %     xlim([-30,30])


end

function print_diffs(x1,x2)

    rmsd_abs=mean_diff('abs',x1,x2,'rmsd');
    rmsd_rel=mean_diff('rel',x1,x2,'rmsd');
    [abs_diff,abs_sigma]=mean_diff('abs',x1,x2);
    [rel_diff,rel_sigma]=mean_diff('rel',x1,x2);

    fprintf('RMSD abs.: %3.3g\n', rmsd_abs)
    fprintf('RMSD rel.: %.1f %%\n', rmsd_rel)
    fprintf('Mean abs. diff.: %3.3g +- %3.3g\n', [abs_diff,abs_sigma])
    fprintf('Mean rel. diff.: %.1f +- %.1f %%\n', [rel_diff,rel_sigma])



end



