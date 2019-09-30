function compare_datasets_diff()
% function to compare various ground-based and satellite datasets and produce
% relative/absolute difference plots
%
% Kristof Bognar, December 2017 

du=2.687e16;

% what to compare
o3_tc=1;
o3_tc_cf=0;
o3_gb=0;
no2_tc=1;
no2_gb=0;

% needs to be run separately from gb comparisons
o3_sat=0;
no2_sat=0;


test=0;

do_violin=0;
violin_type='rel';

% day range for seasonal plots where the end of 'spring' is the last day
% when the sun dips below the horizon
% day_range=[40.25,105.25];
% day_range=[106.25,239.25];
% day_range=[240.25,365.25];
% year_range=[2000,2011];
% year_range=[2000,2006];

% load data
if o3_sat || no2_sat
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
else
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
end

load(data_file)

bruker_no2(bruker_no2.dof_part<1,:)=[];

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

    ind=find(saoz_o3_fixRCD.year<year_range(1)-1 | saoz_o3_fixRCD.year>year_range(2)-1);
    saoz_o3_fixRCD(ind,:)=[];
    
    ind=find(saoz_no2_fixRCD.year<year_range(1)-1 | saoz_no2_fixRCD.year>year_range(2)-1);
    saoz_no2_fixRCD(ind,:)=[];
    
    ind=find(saoz_no2_fixRCD2.year<year_range(1)-1 | saoz_no2_fixRCD2.year>year_range(2)-1);
    saoz_no2_fixRCD2(ind,:)=[];
   
    % bruker and brewer
    ind=find(brewer_o3_ds.year<year_range(1)-1 | brewer_o3_ds.year>year_range(2)-1);
    brewer_o3_ds(ind,:)=[];

    ind=find(bruker_o3.year<year_range(1)-1 | bruker_o3.year>year_range(2)-1);
    bruker_o3(ind,:)=[];

end

%% filter by time
% GBS and saoz
if exist('day_range','var')
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


if ~do_violin
%% difference plots
%% O3 total columns 
if o3_tc
    
    disp('O3 TC')
    
    figure()
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
    disp('PARIS vs OSIRIS O3')
    tmp_time(paris_o3,osiris_o3_bk,label_paris,label_osiris,'');
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

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('PARIS vs ACE-FTS O3 all')
    tmp_time(paris_o3,ace_fts_o3,label_paris,label_fts,'');
    fprintf('\n')

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

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('PARIS vs ACE-MAESTRO O3')
    tmp_time(paris_o3,ace_mae_o3_bw,label_paris,label_mae,'bottom');
    fprintf('\n')
    

    figure()
    disp('Brewer vs OSIRIS O3')
    tmp_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'');
    fprintf('\n')
    
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
    tmp_twilight(ut_o3_cf,ace_mae_o3,'UT_C_F',label_mae,'bottom');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-MAESTRO')
    tmp_twilight(saoz_o3_cf,ace_mae_o3,'SA_C_F',label_mae,'bottom');
    fprintf('\n')

end

%% O3 satellite partial columns
if o3_sat
    
    disp('O3 SAT')
    
    figure(1)
    set(gcf, 'Position', [100, 100, 1200, 650]);    
    
    plotrow=2;
    plotcol=3;
    ii=4;
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
    
    % ACE-FTS vs others
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-MAESTRO vs ACE-FTS')
    tmp_twilight(ace_mae_o3,ace_fts_o3,label_mae,label_fts,'bottomleft', true);
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
    
%     fig_ax = tight_subplot(4,4,[0.05,0.03],[0.1,0.1],[0.1,0.1]);
%     axes(fig_ax(16))
    fig_ax = tight_subplot(3,4,[0.05,0.03],[0.1,0.1],[0.1,0.1]);
    axes(fig_ax(8))
    delete(gca)
    axes(fig_ax(12))
    delete(gca)
    
%     [saoz_o3,saoz_o3_dailyRCD,gbs_o3]=match_saoz_gbs(saoz_o3,saoz_o3_dailyRCD,gbs_o3);
%     [saoz_o3,saoz_o3_dailyRCD2,gbs_o3]=match_saoz_gbs(saoz_o3,saoz_o3_dailyRCD2,gbs_o3);
    
    saoz_o3_fixRCD=saoz_o3_fixRCD2;
    saoz_o3_dailyRCD=saoz_o3_dailyRCD2;
    
%     % brewer on x-axis
%     axes(fig_ax(ii))
%     tmp_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(brewer_o3_ds,saoz_o3_fixRCD,label_brewer,[label_saoz '_{fix}'],'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(brewer_o3_ds,saoz_o3_dailyRCD,label_brewer,[label_saoz '_{day}'],'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'bottomleft',12,false);
%     ii=ii+1;
% 
%     % bruker on x-axis
%     axes(fig_ax(ii))
%     tmp_time(bruker_o3,gbs_o3,label_bruker,label_gbs,'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(bruker_o3,saoz_o3,label_bruker,label_saoz,'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(bruker_o3,saoz_o3_fixRCD,label_bruker,[label_saoz '_{fix}'],'bottomleft',12,false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_time(bruker_o3,saoz_o3_dailyRCD,label_bruker,[label_saoz '_{day}'],'bottomleft',12,false);
%     ii=ii+1;
    
%     % saoz on x-axis
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottomleft',false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3,saoz_o3_fixRCD,label_saoz,[label_saoz '_{fix}'],'bottomleft',false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3,saoz_o3_dailyRCD,label_saoz,[label_saoz '_{day}'],'bottomleft',false);
%     ii=ii+1;
%     
%     % GBS on x-axis
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3_fixRCD,gbs_o3,[label_saoz '_{fix}'],label_gbs,'bottomleft',false);
%     ii=ii+1;
% 
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3_dailyRCD,gbs_o3,[label_saoz '_{day}'],label_gbs,'bottomleft',false);
%     ii=ii+1;
%     
%     % saoz fix vs daily
%     axes(fig_ax(ii))
%     tmp_twilight(saoz_o3_fixRCD,saoz_o3_dailyRCD,[label_saoz '_{fix}'],[label_saoz '_{day}'],'bottomleft',false);
%     ii=ii+1;
    
    % diff plots matching corr plot in paper appendix
    % all minus brewer
    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'',false,1);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'left',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,paris_o3,label_brewer,label_paris,'',false,1);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'',false,2);
    ii=ii+1;
    
    % all minus bruker
    axes(fig_ax(ii))
    tmp_time(bruker_o3,gbs_o3,label_bruker,label_gbs,'',false,1);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(bruker_o3,saoz_o3,label_bruker,label_saoz,'left',false,2);
    ii=ii+1;
    
    axes(fig_ax(ii))
    tmp_time(bruker_o3,paris_o3,label_bruker,label_paris,'',false,11);
    ii=ii+2;
    
    % all minus PARIS
    axes(fig_ax(ii))
    tmp_time(paris_o3,gbs_o3,label_paris,label_gbs,'',false,1);
    ii=ii+1;

    axes(fig_ax(ii))
    tmp_time(paris_o3,saoz_o3,label_paris,label_saoz,'bottomleft',false,2);
    ii=ii+1;
    
    % doas
    axes(fig_ax(ii))
    tmp_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottom',false,11);

end

    
%% NO2
if no2_tc
    
    disp('NO2 TC')
    
    figure()
    plotrow=3;
    plotcol=3;
    ii=1;
    fprintf('\n')

    % OSIRIS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs OSIRIS')
    tmp_time(gbs_no2,osiris_no2,label_gbs,label_osiris,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GU vs OSIRIS')
    tmp_time(gbs_no2uv,osiris_no2uv,label_gbs_uv,label_osiris,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs OSIRIS')
    tmp_time(saoz_no2,osiris_no2,label_saoz,label_osiris,'left');
    fprintf('\n')
    
% % %     % OSIRIS vs SAOZ all year
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('SAOZ fixRCD vs OSIRIS')
% % %     tmp_time(saoz_no2_fixRCD,osiris_no2,'SAOZ_{fix}',label_osiris,'left');
% % %     fprintf('\n')
% % % 
% % %     subplot(plotrow,plotcol,ii)
% % %     ii=ii+1;
% % %     disp('SAOZ dsilyRCD vs OSIRIS')
% % %     tmp_time(saoz_no2_dailyRCD,osiris_no2,'SAOZ_{day}',label_osiris,'left');
% % %     fprintf('\n')

%     ii=ii+1;
    
    % ACE-FTS vs GB
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('UT-GBS vs ACE-FTS all')
    tmp_twilight(gbs_no2,ace_fts_no2,label_gbs,label_fts,'left');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GU vs ACE-FTS all')
    tmp_twilight(gbs_no2uv,ace_fts_no2uv,label_gbs_uv,label_fts,'left');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ vs ACE-FTS all')
    tmp_twilight(saoz_no2,ace_fts_no2,label_saoz,label_fts,'');
    fprintf('\n')

    % Sat vs Bruker
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Bruker vs OSIRIS all')
    tmp_time(bruker_no2,osiris_no2_bk,label_bruker,label_osiris,'');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('Bruker vs ACE-FTS all')
    tmp_time(bruker_no2,ace_fts_no2_bk,label_bruker,label_fts,'');
    fprintf('\n')
    
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

%% NO2 satellite partial columns
if no2_sat
    
    disp('NO2 SAT')
    
    figure()
    plotrow=1;
    plotcol=1;
    ii=1;
    fprintf('\n')
    
%     % OSIRIS vs others
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs OSIRIS')
%     tmp_time(ace_mae_no2,osiris_no2,label_mae,label_osiris,'bottomleft', 12, true);
%     fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('ACE-FTS vs OSIRIS')
    ace_fts_no2.part_col=ace_fts_no2.part_col_os;
    tmp_time(ace_fts_no2,osiris_no2,label_fts,label_osiris,'bottomleft', 12, true);
    fprintf('\n')

    ii=ii+1;
    
%     % ACE-FTS vs others
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs ACE-FTS')
%     tmp_twilight(ace_mae_no2,ace_fts_no2,label_mae,label_fts,'bottomleft', true);
%     fprintf('\n')
% 
%     ii=ii+1;
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-FTS vs OSIRIS')
%     tmp_time(osiris_no2,ace_fts_no2,label_osiris,label_fts,'bottomleft', 12, true);
%     fprintf('\n')
% 
%     % ACE-MAESTRO vs others
%     ii=ii+1;
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs ACE-FTS')
%     tmp_twilight(ace_fts_no2,ace_mae_no2,label_fts,label_mae,'bottomleft', true);
%     fprintf('\n')
%     
%     subplot(plotrow,plotcol,ii)
%     ii=ii+1;
%     disp('ACE-MAESTRO vs OSIRIS')
%     tmp_time(osiris_no2,ace_mae_no2,label_osiris,label_mae,'bottomleft', 12, true);
%     fprintf('\n')
    
end

%% NO2 from ground-based instruments
if no2_gb
    
    disp('NO2 GB')
    
    figure()
    plotrow=2;
    plotcol=3;
    ii=1;
    fprintf('\n')

% %     bruker_no2.tot_col=bruker_no2.part_col_doas;
    bruker_no2.tot_col=bruker_no2.part_col;
%     bruker_no2(2355:end,:)=[];
%     bruker_no2(bruker_no2.sza>80,:)=[];

    
% 4x4 plot using different saoz retrievals
% %     [saoz_no2,saoz_no2_dailyRCD,gbs_no2uv]=match_saoz_gbs(saoz_no2,saoz_no2_dailyRCD,gbs_no2uv);
% %     [saoz_no2,saoz_no2_dailyRCD2,gbs_no2uv]=match_saoz_gbs(saoz_no2,saoz_no2_dailyRCD2,gbs_no2uv);
% % 
% %     [saoz_no2,saoz_no2_dailyRCD,gbs_no2]=match_saoz_gbs(saoz_no2,saoz_no2_dailyRCD,gbs_no2);
% %     [saoz_no2,saoz_no2_dailyRCD2,gbs_no2]=match_saoz_gbs(saoz_no2,saoz_no2_dailyRCD2,gbs_no2);
% % 
% %     saoz_no2_dailyRCD=saoz_no2_dailyRCD2;
% %     saoz_no2_fixRCD=saoz_no2_fixRCD2;
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('UT NO2 vs GBS NO2 UV')
% %     tmp_twilight(gbs_no2,gbs_no2uv,label_gbs,label_gbs_uv,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ NO2 vs GBS NO2 UV')
% %     tmp_twilight(saoz_no2,gbs_no2uv,label_saoz,label_gbs_uv,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ fixRCD vs GBS NO2 UV')
% %     tmp_twilight(saoz_no2_fixRCD,gbs_no2uv,'SAOZ_{fix}',label_gbs_uv,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ dailyRCD vs GBS NO2 UV')
% %     tmp_twilight(saoz_no2_dailyRCD,gbs_no2uv,'SAOZ_{day}',label_gbs_uv,'bottomleft');
% %     fprintf('\n')
% %     
% %     %%%
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('BK vs GBS NO2 UV')
% %     tmp_time(bruker_no2,gbs_no2uv,label_bruker,label_gbs_uv,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ NO2 vs GBS NO2')
% %     tmp_twilight(saoz_no2,gbs_no2,label_saoz,label_gbs,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ fixRCD vs GBS NO2')
% %     tmp_twilight(saoz_no2_fixRCD,gbs_no2,'SAOZ_{fix}',label_gbs,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ dailyRCD vs GBS NO2')
% %     tmp_twilight(saoz_no2_dailyRCD,gbs_no2,'SAOZ_{day}',label_gbs,'bottomleft');
% %     fprintf('\n')
% %     
% %     %%%
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('BK vs GBS NO2')
% %     tmp_time(bruker_no2,gbs_no2,label_bruker,label_gbs,'bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ fixRCD vs SAOZ dailyRCD')
% %     tmp_twilight(saoz_no2_fixRCD,saoz_no2_dailyRCD,'SAOZ_{fix}','SAOZ_{day}','bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ vs SAOZ fixRCD')
% %     tmp_twilight(saoz_no2,saoz_no2_fixRCD,label_saoz,'SAOZ_{fix}','bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('SAOZ vs SAOZ dailyRCD')
% %     tmp_twilight(saoz_no2,saoz_no2_dailyRCD,label_saoz,'SAOZ_{day}','bottomleft');
% %     fprintf('\n')
% % 
% %     %%%
% %     
% %     ii=ii+1;
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('BK vs SAOZ')
% %     tmp_time(bruker_no2,saoz_no2,label_bruker,label_saoz,'bottomleft');
% %     fprintf('\n')
% %     
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('BK vs SAOZ fixRCD')
% %     tmp_time(bruker_no2,saoz_no2_fixRCD,label_bruker,'SAOZ_{fix}','bottomleft');
% %     fprintf('\n')
% % 
% %     subplot(plotrow,plotcol,ii)
% %     ii=ii+1;
% %     disp('BK vs SAOZ dailyRCD')
% %     tmp_time(bruker_no2,saoz_no2_dailyRCD,label_bruker,'SAOZ_{day}','bottomleft');
% %     fprintf('\n')


% to match corr plot in paper appendix

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ NO2 vs GBS NO2')
    tmp_twilight(saoz_no2,gbs_no2,label_saoz,label_gbs,'bottomleft');
    fprintf('\n')

    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GBS NO2 vis vs Bruker NO2')
    tmp_time(gbs_no2,bruker_no2,label_gbs,label_bruker,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GBS NO2 vis vs GBS NO2 UV')
    tmp_twilight(gbs_no2,gbs_no2uv,label_gbs,label_gbs_uv,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ NO2 vs GBS NO2 UV')
    tmp_twilight(saoz_no2,gbs_no2uv,label_saoz,label_gbs_uv,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('GBS NO2 UV vs Bruker NO2')
    tmp_time(gbs_no2uv,bruker_no2,label_gbs_uv,label_bruker,'bottomleft');
    fprintf('\n')
    
    subplot(plotrow,plotcol,ii)
    ii=ii+1;
    disp('SAOZ NO2 vs Bruker NO2')
    tmp_time(saoz_no2,bruker_no2,label_saoz,label_bruker,'bottomleft');
    fprintf('\n')
    
end

if test
    

figure();

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


else % do violin plots
    
label_size=15;
    
if o3_tc    
        
    figure()

    subplot(211)
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_o3,  gbs_o3,  false );
    af_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_o3, saoz_o3,  false );
    af_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(ace_fts_o3_bk, bruker_o3, 12,  false );
    af_bk=rel_abs_diff(violin_type,x1,x2); 

    [ x1, x2 ] = find_coincidences_twilight(ace_mae_o3,  gbs_o3,  false );
    am_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_mae_o3, saoz_o3,  false );
    am_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(ace_mae_o3_bk, bruker_o3, 12,  false );
    am_bk=rel_abs_diff(violin_type,x1,x2); 
 
    
    diff={af_gv,af_sa,af_bk,am_gv,am_sa,am_bk};

    % cristen's values
    if strcmp(violin_type,'abs')
        c_diffs=[28.4,19.2,-21.6,21.7,8.7,-29.8];
    else
        c_diffs=[6.5,4.8,-4.7,5,1.6,-6.1];
    end
    c_diffs=[1,1,1,1,1,1]*1000;
    
    violin_mod(diff,'mc','k','medc',c_diffs,...
               'edgecolor','','facecolor',[0 0 0.2],'facealpha',0.3,...
               'xlabel',{'FTS-GBS','FTS-SAOZ','FTS-Bruker','MAE-GBS','MAE-SAOZ','MAE-Bruker'});
%     set(gcf, 'Position', [100, 100, 1200, 700]);    

%%% repeat with smoothing

    load('/home/kristof/work/satellite_validation/all_data_smooth.mat');
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_o3,  gbs_o3,  false );
    af_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_o3, saoz_o3,  false );
    af_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(ace_fts_o3_bk, bruker_o3, 12,  false );
    af_bk=rel_abs_diff(violin_type,x1,x2); 

    [ x1, x2 ] = find_coincidences_twilight(ace_mae_o3,  gbs_o3,  false );
    am_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_mae_o3, saoz_o3,  false );
    am_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(ace_mae_o3_bk, bruker_o3, 12,  false );
    am_bk=rel_abs_diff(violin_type,x1,x2); 
 
    diff={af_gv,af_sa,af_bk,am_gv,am_sa,am_bk};
    
    violin_mod2(diff,'mc','r--','medc',c_diffs,...
               'edgecolor','r','facecolor',[1 0 0.2],'facealpha',0,...
               'xlabel',{'FTS-GBS','FTS-SAOZ','FTS-Bruker','MAE-GBS','MAE-SAOZ','MAE-Bruker'});

%%%           
           
           
    if strcmp(violin_type,'abs')       
        ylim([-150,150])
        ylabel('\Delta_{abs} (DU)')
    else
        ylim([-50,50])
        ylabel('\Delta_{rel} (%)')
    end
    
    ax = gca;
    ax.YGrid = 'on';
    set(ax,'yminorgrid','on')
 
%%%%%%%%%%%%%%%%%    
    
    subplot(212)
    load('/home/kristof/work/satellite_validation/all_data_nosmooth.mat');
    [ x1, x2 ] = find_coincidences_time(osiris_o3,  gbs_o3, 12, false );
    os_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_time(osiris_o3, saoz_o3, 12, false );
    os_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(osiris_o3_bk, bruker_o3, 12,  false );
    os_bk=rel_abs_diff(violin_type,x1,x2); 

    [ x1, x2 ] = find_coincidences_time(osiris_o3_bw,  brewer_o3_ds, 12, false );
    os_bw=rel_abs_diff(violin_type,x1,x2);
    
    diff={os_gv,os_sa,os_bk,os_bw};

    if strcmp(violin_type,'abs')
        c_diffs=[19.9,32.2,-0.4,10.3];
    else
        c_diffs=[5.7,7.3,0.1,2.8];
    end
    c_diffs=[1,1,1,1]*1000;
    
    violin_mod(diff,'mc','k','medc',c_diffs,...
               'edgecolor','','facecolor',[0 0 0.2],'facealpha',0.3,...
               'xlabel',{'OS-GBS','OS-SAOZ','OS-Bruker','OS-Brewer'});
    set(gcf, 'Position', [100, 100, 1000, 700]);    

%%% repeat for smoothed

    load('/home/kristof/work/satellite_validation/all_data_smooth.mat');
    [ x1, x2 ] = find_coincidences_time(osiris_o3,  gbs_o3, 12, false );
    os_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_time(osiris_o3, saoz_o3, 12, false );
    os_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(osiris_o3_bk, bruker_o3, 12,  false );
    os_bk=rel_abs_diff(violin_type,x1,x2); 

    [ x1, x2 ] = find_coincidences_time(osiris_o3_bw,  brewer_o3_ds, 12, false );
    os_bw=rel_abs_diff(violin_type,x1,x2);
    
    diff={os_gv,os_sa,os_bk,os_bw};

    violin_mod2(diff,'mc','r--','medc',c_diffs,...
               'edgecolor','r','facecolor',[1 0 0.2],'facealpha',0,...
               'xlabel',{'OS-GBS','OS-SAOZ','OS-Bruker','OS-Brewer'});

%%%    
    
    if strcmp(violin_type,'abs')       
        ylim([-150,150])
        ylabel('\Delta_{abs} (DU)')
    else
        ylim([-50,50])
        ylabel('\Delta_{rel} (%)')
    end
    
    ax = gca;
    ax.YGrid = 'on';
    set(ax,'yminorgrid','on')

end


if no2_tc
    
    figure()

    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2,  gbs_no2,  false );
    af_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2uv, gbs_no2uv,  false );
    af_gu=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2, saoz_no2,  false );
    af_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(osiris_no2,  gbs_no2, 12, false );
    os_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_time(osiris_no2uv,  gbs_no2uv, 12, false );
    os_gu=rel_abs_diff(violin_type,x1,x2);

    [ x1, x2 ] = find_coincidences_time(osiris_no2, saoz_no2, 12, false );
    os_sa=rel_abs_diff(violin_type,x1,x2);    
    

    diff={af_gv,af_gu,af_sa,os_gv,os_gu,os_sa};

    if strcmp(violin_type,'abs')
%         c_diffs=[1.6,1.7,1.7,-3,-2,2]*1e15;
        c_diffs=[1.6,1.7,1.7,-3,-2,2]*1e14;
    else
        c_diffs=[15.2,13.6,12.7,-7.3,-3.3,10.2];
    end
    c_diffs=[1,1,1,1,1,1]*1000;
   
    violin_mod(diff,'mc','k','medc',c_diffs,...
               'edgecolor','','facecolor',[0 0 0.2],'facealpha',0.3,...
               'xlabel',{'FTS-GBS','FTS-GBS UV ','FTS-SAOZ','OS-GBS','OS-GBS UV','OS-SAOZ'});

%%% repeat for smoothed
    load('/home/kristof/work/satellite_validation/all_data_smooth.mat');

    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2,  gbs_no2,  false );
    af_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2uv, gbs_no2uv,  false );
    af_gu=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_twilight(ace_fts_no2, saoz_no2,  false );
    af_sa=rel_abs_diff(violin_type,x1,x2);    
    
    [ x1, x2 ] = find_coincidences_time(osiris_no2,  gbs_no2, 12, false );
    os_gv=rel_abs_diff(violin_type,x1,x2);
    
    [ x1, x2 ] = find_coincidences_time(osiris_no2uv,  gbs_no2uv, 12, false );
    os_gu=rel_abs_diff(violin_type,x1,x2);

    [ x1, x2 ] = find_coincidences_time(osiris_no2, saoz_no2, 12, false );
    os_sa=rel_abs_diff(violin_type,x1,x2);    
    

    diff={af_gv,af_gu,af_sa,os_gv,os_gu,os_sa};

    violin_mod2(diff,'mc','r--','medc',c_diffs,...
               'edgecolor','r','facecolor',[1 0 0.2],'facealpha',0,...
               'xlabel',{'FTS-GBS','FTS-GBS UV ','FTS-SAOZ','OS-GBS','OS-GBS UV','OS-SAOZ'});

%%%           
           
    if strcmp(violin_type,'abs')
        ylim([-3.1,3.1]*1e15)
        ylabel('\Delta_{abs} (molec/cm^2)')
    else
        ylim([-100,100])
        ylabel('\Delta_{rel} (%)')
    end
    
    set(gcf, 'Position', [100, 100, 1200, 700]);

    
end
set(findall(gcf,'-property','FontSize'),'FontSize',label_size)

ax = gca;
ax.YGrid = 'on';
set(ax,'yminorgrid','on')

end

end

%% local functions to avoid copy-pasting code
function tmp_twilight(table1, table2, instr1, instr2, loc_str, partcol,xx)

    if nargin<6, partcol=false; end
    
    [ x1, x2, e1, e2, times, DMP] =...
        find_coincidences_twilight( table2, table1, partcol );
    
%     T1=DMP.T1(:,14);
%     T2=DMP.T2(:,14);
%     spv2=DMP.spv2(:,14);
    
    sza1=DMP.sza1; % sat tangent point
    sza2=DMP.sza2; % g-b SZA
    
%     diff_plot_time('abs', x1,x2, times, 'timeseries', instr2, instr1)
%     diff_plot_time('abs', x1,x2, times, 'seasonal', instr2, instr1)
%     diff_plot_time('rel', x1,x2, times, 'ymean', instr2, instr1)
     diff_plot_param('abs', x1,x2, sza2, 'SZA', instr2, instr1)
%      diff_plot_param('rel', x1,x2, T1, 'T', instr2, instr1)
%      diff_plot_param('rel', x1,x2, spv2, 'sPV', instr2, instr1)
%      diff_plot_param('abs', x1,x2, (x1+x2)/2, 'mean', instr2, instr1)
%      diff_plot_param('abs', x1,x2, mjd2k_to_ft(times), 'day', instr2, instr1)

    abs_diff=mean_diff('abs',x1,x2);
    rel_diff=mean_diff('rel',x1,x2);

    fprintf('Mean abs. diff.: %3.3g\n', abs_diff)
    fprintf('Mean rel. diff.: %.1f %%\n', rel_diff)

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

function tmp_time(table1, table2, instr1, instr2, loc_str, partcol, xx)

    if nargin==5, partcol=false; end
    dt=12;
    
    [ x1, x2, e1, e2, times,~, DMP ] =...
        find_coincidences_time( table2, table1, dt, partcol );
  
%     T1=DMP.T1(:,14);
%     T2=DMP.T2(:,14);
%     spv2=DMP.spv2(:,14);
    
    sza1=DMP.sza1; % sat tangent point
    sza2=DMP.sza2; % g-b SZA
        
%     diff_plot_time('abs', x1,x2, times, 'timeseries', instr2, instr1)
%     diff_plot_time('abs', x1,x2, times, 'seasonal', instr2, instr1)
%     diff_plot_time('rel', x1,x2, times, 'ymean', instr2, instr1)
    diff_plot_param('abs', x1,x2, sza2, 'SZA', instr2, instr1)
%     diff_plot_param('rel', x1,x2, T2, 'T', instr2, instr1)
%     diff_plot_param('rel', x1,x2, spv2, 'sPV', instr2, instr1)
%     diff_plot_param('abs', x1,x2, (x1+x2)/2, 'mean', instr2, instr1)
%     diff_plot_param('abs', x1,x2, mjd2k_to_ft(times), 'day', instr2, instr1)
     
    abs_diff=mean_diff('abs',x1,x2);
    rel_diff=mean_diff('rel',x1,x2);

    fprintf('Mean abs. diff.: %3.3g\n', abs_diff)
    fprintf('Mean rel. diff.: %.1f  %%\n', rel_diff)

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

function [saoz,saoz_dailyRCD,gbs]=match_saoz_gbs(saoz,saoz_dailyRCD,gbs)

    % match saoz V3 data to saoz daily RCD
    arr1=saoz(:,1:3);
    arr2=saoz_dailyRCD(:,1:3);
    [~,ind1,ind2]=intersect(arr1,arr2);
    saoz=saoz(ind1,:);
    saoz_dailyRCD=saoz_dailyRCD(ind2,:);
    
    % match to GBS
    arr1=saoz(:,1:3);
    arr2=gbs(:,1:3);
    [~,ind1,ind2]=intersect(arr1,arr2);
    saoz=saoz(ind1,:);
    saoz_dailyRCD=saoz_dailyRCD(ind1,:);
    gbs=gbs(ind2,:);
    
    
end
