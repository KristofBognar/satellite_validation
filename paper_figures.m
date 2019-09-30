function paper_figures()
%PAPER_FIGURES plot all figures for sat val paper

data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';

% data_file='/home/kristof/work/satellite_validation/all_data_nosmooth_BK_full_corr.mat';
% ind=find(bruker_o3.year==2014 | bruker_o3.year==2015);
% bruker_o3=bruker_o3(ind,:);

load(data_file)

%% select figures/settings

% save figures?
save_figs=1;
%

meas_loc=0;
sza_avk=0;

diff_hist=0;

o3_sat_sat=0;

o3_OS_BW=0;
o3_corr=0;
o3_diff=0;

o3_prof_diff=0;

no2_corr=0;
no2_diff=1;

smooth_diffs=0;

tc_stuff=0;

o3_gb=0;
no2_gb=0;


%% global variables
global alphabet text_size label_size ii alt_dmp_in dmp_filter save_as_pdf

text_size=9; % matlab default is 8
label_size=10.5; % matlab default is 10
alphabet = 'abcdefghijklmnopqrstuvwxyz';
alt_dmp_in=alt_dmp;

save_as_pdf=1;

% use tighter sprintime coincidence criteria?
% 1 for sPV and T limits, 2 for latitude limit
dmp_filter=0;

% partial coulon comparison with Bruker?
bruker_part_col=0;

% filter datasets to match Cristen's dates
cristen=0;

% exclude summer data
no_summer=0;

% filter by season
pick_season=0;

% extra filters
% % gbs_o3(gbs_o3.year<2010,:)=[];
% % saoz_o3(saoz_o3.year<2010,:)=[];
% % brewer_o3_ds(brewer_o3_ds.sza>76.35,:)=[];

% kick out lamp signal data for UT-GBS in 2017
% gbs_o3(3005:3117,:)=[];
% gbs_no2(3012:3123,:)=[];

% kick out bad shutter data for UT-GBS in 2017
% gbs_o3(gbs_o3.year==2017 & gbs_o3.fractional_time<71,:)=[];
% gbs_no2(gbs_no2.year==2017 & gbs_no2.fractional_time<71,:)=[];

bruker_no2.tot_col=bruker_no2.part_col; % use 12-40 km col for comparisons
% bruker_no2(bruker_no2.sza>80,:)=[];
% bruker_no2(bruker_no2.dof_part<1,:)=[];

if pick_season
    
    % day ranges, end days are included in range
    day_range=[40.25,105.25]; disp('Season restricted to spring');
%     day_range=[106.25,239.25]; disp('Season restricted to summer');
%     day_range=[240.25,365.25]; disp('Season restricted to fall');

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
    
    ind=find(bruker_no2.fractional_time<day_range(1)-1 | bruker_no2.fractional_time>day_range(2));
    bruker_no2(ind,:)=[];
    
    ind=find(paris_o3.fractional_time<day_range(1)-1 | paris_o3.fractional_time>day_range(2));
    paris_o3(ind,:)=[];

end

if no_summer
    
    day_range=[106.25,239.25]; disp('Summer data excluded');
%     day_range=[80.25,259.25]; disp('Summer data excluded');
    
    ind=find(gbs_o3.fractional_time>day_range(1)-1 & gbs_o3.fractional_time<day_range(2));
    gbs_o3(ind,:)=[];

    ind=find(gbs_no2.fractional_time>day_range(1)-1 & gbs_no2.fractional_time<day_range(2));
    gbs_no2(ind,:)=[];

    ind=find(gbs_no2uv.fractional_time>day_range(1)-1 & gbs_no2uv.fractional_time<day_range(2));
    gbs_no2uv(ind,:)=[];

    ind=find(saoz_o3.fractional_time>day_range(1)-1 & saoz_o3.fractional_time<day_range(2));
    saoz_o3(ind,:)=[];

    ind=find(saoz_no2.fractional_time>day_range(1)-1 & saoz_no2.fractional_time<day_range(2));
    saoz_no2(ind,:)=[];

    ind=find(saoz_no2_fixRCD.fractional_time>day_range(1)-1 &...
             saoz_no2_fixRCD.fractional_time<day_range(2));
    saoz_no2_fixRCD(ind,:)=[];

    % bruker and brewer
    ind=find(brewer_o3_ds.fractional_time>day_range(1)-1 & brewer_o3_ds.fractional_time<day_range(2));
    brewer_o3_ds(ind,:)=[];

    ind=find(bruker_o3.fractional_time>day_range(1)-1 & bruker_o3.fractional_time<day_range(2));
    bruker_o3(ind,:)=[];

    ind=find(bruker_no2.fractional_time>day_range(1)-1 & bruker_no2.fractional_time<day_range(2));
    bruker_no2(ind,:)=[];
    
end


if cristen % filter data based on times in Cristen's paper for direct comparisons
    
    disp('Dates matched to Adams et al. (2012)');
    
    gbs_o3(1442:end,:)=[];
    saoz_o3(567:end,:)=[];
    brewer_o3_ds(28772:end,:)=[];
    bruker_o3(2517:end,:)=[];
    ace_mae_o3(185:end,:)=[];
    ace_fts_o3(289:end,:)=[];
    
    gbs_no2(1420:end,:)=[];
    gbs_no2uv(739:end,:)=[];
    saoz_no2(567:end,:)=[];
    bruker_no2(2355:end,:)=[];
    osiris_no2(1775:end,:)=[];
    osiris_no2uv(1775:end,:)=[];
    
end

if bruker_part_col, disp('Bruker/PARIS partial columns used'); end

if dmp_filter==1
    disp('DMP filter used for coincidences');
elseif dmp_filter==2
    disp('Lat filter used for coincidences');
end


%% Measurement locations on map
if meas_loc

    plot_geoloc;

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 390]);
    
    save_pdf(save_figs,'geoloc')


end

%% histogram of relative differences
if diff_hist

%     plot_hist_adams_fraser_griffin(1);
% 
%     set(findall(gcf,'-property','FontName'),'FontName','Arial')
% %     set(gca, 'FontSize', label_size)
%     set(gca, 'FontSize', 14)
%     set(gcf, 'Position', [100, 100, 1000, 600]);
%     
%     save_pdf(save_figs,'hist_o3')
% 
%     plot_hist_adams_fraser_griffin(2);
% 
%     set(findall(gcf,'-property','FontName'),'FontName','Arial')
% %     set(gca, 'FontSize', label_size)
%     set(gca, 'FontSize', 14)
%     set(gcf, 'Position', [100, 100, 800, 600]);
% 
%     save_pdf(save_figs,'hist_no2')

    plot_hist_adams_fraser_griffin(99);

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
%     set(gca, 'FontSize', label_size+2)
%     set(gca, 'FontSize', 14)
    set(gcf, 'Position', [100, 100, 1000, 750]);
    
    save_pdf(save_figs,'hist_all')
    
end
%% yearly SZA variation and spring/summer AVKs
if sza_avk
    
    figure()
    
    ax1=subplot(3,3,[1 2 3]);
    
    t1 = datetime(2017,2,11,0,0,0);
    t2 = datetime(2017,10,29,0,0,0);
    date_tmp=[t1:minutes(30):t2]';
    
    [~, sza] = SolarAzEl(date_tmp,ones(size(date_tmp))*80.05,ones(size(date_tmp))*-86.42,...
                             0.61*ones(size(date_tmp)));
    sza=90-sza;    
    date_tmp(1:50)=date_tmp(50);
    date_tmp(end-49:end)=date_tmp(end-49);
    
    patch([0 1 1 0],[0 0 1 1],[0 .6 .6],'edgecolor','none'), hold on
    plot([0,1],[91,91],'k-')
    
    plot(date_tmp,sza,'.','color',[0 .6 .6]), hold on
    plot([date_tmp(1),date_tmp(end)],[91,91],'k-')
    plot([date_tmp(1),date_tmp(end)],[86,86],'k-')

    xlim(datenum([datetime(2017,2,10,19,0,0),datetime(2017,10,29,9,0,0)]))
    ylim([50,120])
    ylabel('SZA')
    
    ll=legend('Available SZA','Ideal DOAS SZA range','orientation','horizontal');
    set(ll,'position',[0.26 0.894 0.5, 0.04],'box','off')
    
    datetick('x','mmm','keeplimits','keepticks')
    
    
    % get DOAS AVKs
    year_tmp=[2017;2017];
    % spring equinox, SZA=90 and summer solstice, SZA=76.5 in 2017 (March 20, June 21)
    ft_tmp=[78.494;171.244];
    
    AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';
    load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
    % time when we have the sonde measurement
    sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
    sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)
    % sat measurement time
    measurement_time = year_tmp + ft_tmp./daysinyear(year_tmp);
    sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

    [~,o3_avk]=read_DOAS_prof_avk(1,[year_tmp,floor(ft_tmp)+1,sonde_in],AVK_LUT_dir );
    [~,no2_avk]=read_DOAS_prof_avk(2,[year_tmp,ft_tmp],AVK_LUT_dir );
    [~,no2uv_avk]=read_DOAS_prof_avk(3,[year_tmp,ft_tmp],AVK_LUT_dir );

    % bruker: March 21, 2017, 16:27 and June 21, 2017, 18:48
    [bk_alt_o3,~,~,bk_avk_o3]=read_bruker_prof_avk(1,[6.2896677199e3,6.3818044792e3]);
    % times are withing 30 min for NO2
    [bk_alt_no2,~,~,bk_avk_no2]=read_bruker_prof_avk(2,[6.2896677199e3,6.3818044792e3]);
    
    alt_km=0.5:59.5;
    
    text(0.05,0.95,'a)','Units','normalized')
    
    box off

    ax21=subplot(334);
    plot(o3_avk(1,4:end),alt_km,'-','color',[0 .6 .6],'linewidth',2), hold on
    plot(no2_avk(1,3:end),alt_km,'-','color',[.6 0 .6],'linewidth',2)
    plot(no2uv_avk(1,3:end),alt_km,'-','color',[.6 .6 0],'linewidth',2)
    plot([1,1],[0,60],'k--')
    xlim([0,1.6])
    ylim([0,60])
    
    text(0.05,0.95,'b)','Units','normalized')
    text(0.05,0.5,'Mar 20','Units','normalized')
    ylabel('Altitude (km)')
    box off    
    
    ax22=subplot(335);
    plot(o3_avk(2,4:end),alt_km,'-','color',[0 .6 .6],'linewidth',2), hold on
    plot(no2_avk(2,3:end),alt_km,'-','color',[.6 0 .6],'linewidth',2)
    plot(no2uv_avk(2,3:end),alt_km,'-','color',[.6 .6 0],'linewidth',2)
    plot([1,1],[0,60],'k--')
    xlim([0,1.6])
    ylim([0,60])
    
    text(0.05,0.5,'Jun 21','Units','normalized')
    text(1.13,0.9,'DOAS AVK','Units','normalized')
    ll=legend('O_3','NO_2-vis','NO_2-UV');
    set(ll,'position',[0.65 0.4 0.25 0.2],'box','off')
    box off    
    

    ax31=subplot(337);
    plot(bk_avk_o3(1,:),bk_alt_o3,'-','color',[0 .6 .6],'linewidth',2), hold on
    plot(bk_avk_no2(1,:),bk_alt_no2,'-','color',[.6 0 .6],'linewidth',2), hold on
    plot([1,1],[0,60],'k--')
    xlim([0,1.6])
    ylim([0,60])
    
    text(0.05,0.95,'c)','Units','normalized')
    text(0.05,0.5,'Mar 20','Units','normalized')
    ylabel('Altitude (km)')
    xlabel('Column AVK')
    box off    
    
    ax32=subplot(338);
    plot(bk_avk_o3(2,:),bk_alt_o3,'-','color',[0 .6 .6],'linewidth',2), hold on
    plot(bk_avk_no2(2,:),bk_alt_no2,'-','color',[.6 0 .6],'linewidth',2), hold on
    plot([1,1],[0,60],'k--')
    xlim([0,1.6])
    ylim([0,60])
    
    text(0.05,0.5,'Jun 21','Units','normalized')
    text(1.13,0.9,'Bruker AVK','Units','normalized')
    xlabel('Column AVK')
    ll=legend('O_3','NO_2');
    set(ll,'position',[0.65 0.11 0.2 0.2],'box','off')
    
    box off    
    
    %%%
    
    set(ax22,'YTickLabel',[]);    
    set(ax32,'YTickLabel',[]);    

    set(ax21,'XTickLabel',[]);    
    set(ax22,'XTickLabel',[]);    
%     set(ax32,'XTickLabel',[{''},{'0.5'},{'1'},{'1.5'}]);    
    set(ax32,'XTickLabel',{'','0.5','1','1.5'});    
    
    set(ax21,'Position',[0.13 0.4 0.25 0.24])
    set(ax22,'Position',[0.38 0.4 0.25 0.24])

    set(ax31,'Position',[0.13 0.11 0.25 0.24])
    set(ax32,'Position',[0.38 0.11 0.25 0.24])
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 500]);
    
    save_pdf(save_figs,'sza_avk')

end

%% sat-sat ozone 
if o3_sat_sat
    
    figure()
    fig_ax = tight_subplot(2,3,[0.18,0.05],[0.11,0.07],[0.11,0.07]);
    ii=1;
    
    % OSIRIS vs ACE-FTS correlation
    axes(fig_ax(ii))
    corr_time(ace_fts_o3,osiris_o3,label_fts,label_osiris,'bottomleft',true);
    ii=ii+1;


    % OSIRIS vs ACE-MAESTRO correlation
    axes(fig_ax(ii))
    corr_time(ace_mae_o3,osiris_o3,label_mae,label_osiris,'bottomleft',true);
    ii=ii+1;
    
    % ACE-FTS vs ACE-MAESTRO correlation
    axes(fig_ax(ii))
    corr_twilight(ace_mae_o3,ace_fts_o3,label_mae,label_fts,'bottomleft',true);
    ii=ii+1;
    
    % OSIRIS vs ACE-FTS differences
    axes(fig_ax(ii))
    diff_time(ace_fts_o3,osiris_o3,label_fts,label_osiris,'bottomleft',true,9);
    ii=ii+1;

    % OSIRIS vs ACE-MAESTRO differences
    axes(fig_ax(ii))
    diff_time(ace_mae_o3,osiris_o3,label_mae,label_osiris,'bottom',true,9);
    ii=ii+1;
    
    % ACE-FTS vs ACE-MAESTRO differences
    axes(fig_ax(ii))
    diff_twilight(ace_mae_o3,ace_fts_o3,label_mae,label_fts,'bottom',true,9);
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 530]);
    
    save_pdf(save_figs,'sat_sat_o3')

end

%% OSIRIS vs brewer ozone columns, correlation and difference plots
if o3_OS_BW
    
    figure()
    
    fig_ax = tight_subplot(2,1,[0.18,0.05],[0.11,0.07],[0.18,0.10]);
    
    ii=1;
    
    axes(fig_ax(ii))
    corr_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'bottomleft',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'bottomleft',false,12);

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 530]);
    
    save_pdf(save_figs,'OS_BW_o3')
    
end

if 0 % OS-BW SZA plot, single figure since it soesn't fit on main ozone diffs plot
    
    figure()
    fig_ax = tight_subplot(1,1,[0.18,0.05],[0.22,0.18],[0.18,0.10]);
    ii=1;
    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'bottomleft',false,12);

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 280]);
    xlim([55,85])
    
    save_pdf(save_figs,'OS_BW_o3_sza')
    
end

%% ozone gb vs sat total column correlation plots
if o3_corr
    
    figure()
    
    fig_ax = tight_subplot(3,4,[0.065,0.03],[0.08,0.06],[0.09,0.05]);

%     axes(fig_ax(8))
%     delete(gca)
%     axes(fig_ax(12))
%     delete(gca)

    ii=1;
    
    % OSIRIS vs GB
    axes(fig_ax(ii))
    corr_time(gbs_o3,osiris_o3,label_gbs,label_osiris,'left',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(saoz_o3,osiris_o3,label_saoz,label_osiris,'',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_time(bruker_o3,osiris_o3_bk,label_bruker,label_osiris,'',bruker_part_col);
    ii=ii+1;

%     axes(fig_ax(ii))
%     corr_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'bottom',false);
%     ii=ii+1;
    axes(fig_ax(ii))
    corr_time(paris_o3,osiris_o3_bk,label_paris,label_osiris,'',bruker_part_col);
    ii=ii+1;

    % ACE-FTS vs GB
    axes(fig_ax(ii))
    corr_twilight(gbs_o3,ace_fts_o3,label_gbs,label_fts,'left',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_twilight(saoz_o3,ace_fts_o3,label_saoz,label_fts,'',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(bruker_o3,ace_fts_o3_bk,label_bruker,label_fts,'',bruker_part_col);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(paris_o3,ace_fts_o3_bk,label_paris,label_fts,'',bruker_part_col);
    ii=ii+1;


    % ACE-MAESTRO vs GB
    axes(fig_ax(ii))
    corr_twilight(gbs_o3,ace_mae_o3,label_gbs,label_mae,'bottomleft',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_twilight(saoz_o3,ace_mae_o3,label_saoz,label_mae,'bottom',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(bruker_o3,ace_mae_o3_bk,label_bruker,label_mae,'bottom',bruker_part_col);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(paris_o3,ace_mae_o3_bk,label_paris,label_mae,'bottom',bruker_part_col);
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 750]);

    save_pdf(save_figs,'corr_o3')

end

%% ozone gb vs sat total column difference plots
if o3_diff
    
    figure()
    
    if save_as_pdf
        fig_ax = tight_subplot(3,2,[0.11,0.05],[0.04,0.05],[0.1,0.05]);
    else
        fig_ax = tight_subplot(3,2,[0.11,0.05],[0.06,0.07],[0.1,0.05]);
    end

    ii=1;
    
    % OSIRIS vs GB
    axes(fig_ax(ii))
    diff_time(gbs_o3,osiris_o3,label_gbs,label_osiris,'',false,1);

    axes(fig_ax(ii))
    diff_time(saoz_o3,osiris_o3,label_saoz,label_osiris,'left',false,2);
    ii=ii+1;
    
    axes(fig_ax(ii))
    diff_time(paris_o3,osiris_o3_bk,label_paris,label_osiris,'',bruker_part_col,1);

    axes(fig_ax(ii))
    diff_time(bruker_o3,osiris_o3_bk,label_bruker,label_osiris,'',bruker_part_col,2);
    ii=ii+1;
    
%     axes(fig_ax(ii))
%     diff_time(brewer_o3_ds,osiris_o3_bw,label_brewer,label_osiris,'',false,2);
%     ii=ii+1;

    % ACE-FTS vs GB
    axes(fig_ax(ii))
    diff_twilight(gbs_o3,ace_fts_o3,label_gbs,label_fts,'',false,1);

    axes(fig_ax(ii))
    diff_twilight(saoz_o3,ace_fts_o3,label_saoz,label_fts,'left',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    diff_time(paris_o3,ace_fts_o3_bk,label_paris,label_fts,'',bruker_part_col, 1);
    
    axes(fig_ax(ii))
%     diff_time(bruker_o3,ace_fts_o3_bk,label_bruker,label_fts,'',bruker_part_col, 11);
    diff_time(bruker_o3,ace_fts_o3_bk,label_bruker,label_fts,'',bruker_part_col, 2);
    ii=ii+1;
    
    % ACE-MAESTRO vs GB
    axes(fig_ax(ii))
    diff_twilight(gbs_o3,ace_mae_o3,label_gbs,label_mae,'',false,1);
    
    axes(fig_ax(ii))
    diff_twilight(saoz_o3,ace_mae_o3,label_saoz,label_mae,'bottomleft',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    diff_time(paris_o3,ace_mae_o3_bk,label_paris,label_mae,'bottom',bruker_part_col,1);

    axes(fig_ax(ii))
    diff_time(bruker_o3,ace_mae_o3_bk,label_bruker,label_mae,'bottom',bruker_part_col,2);
    

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 850]);
    
    save_pdf(save_figs,'diff_o3')

    
end

%% ozone profile difference plots
if o3_prof_diff
    
    
%     % additional filter for MAESTRO based on reprted errors
%     err=(ace_mae_o3.part_col_err./ace_mae_o3.part_col)*100;
%     ace_mae_o3(err>10,:)=[];
%     ace_mae_o3_bk(err>10,:)=[];
%     ace_mae_o3_bw(err>10,:)=[];
    
    % sat-sat profiles
    figure()
    
    if save_as_pdf
        fig_ax = tight_subplot(3,3,[0.07,0.0],[0.08,0.02],[0.13,0.07]);
    else
        fig_ax = tight_subplot(3,3,[0.07,0.0],[0.10,0.04],[0.13,0.07]);
    end

    ii=1;
    diff_prof_time(fig_ax(1:3),osiris_o3,ace_fts_o3,label_osiris,label_fts,...
                       'left', 0)

    ii=2;
    diff_prof_time(fig_ax(4:6),osiris_o3,ace_mae_o3,label_osiris,label_mae,...
                       'left', 0)

    ii=3;
    diff_prof_twilight(fig_ax(7:9),ace_fts_o3,ace_mae_o3,label_fts,label_mae,...
                       'bottomleft', true, 0)
                   
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 850]);
    
    save_pdf(save_figs,'diff_prof_sat_o3')

    % sat vs Bruker profiles -- use Bruker altitude grid
    % smooth with column avk if required
    
    bruker_o3.num_dens=[];
    bruker_o3.Properties.VariableNames{'num_dens_native'} = 'num_dens';
    
    [alt_bk,~,apriori,avk,part_prof,~]=read_bruker_prof_avk(1,osiris_o3_bk.mjd2k,'bruker');
    osiris_o3_bk.num_dens2=NaN(length(osiris_o3_bk.tot_col),20);
    for i=1:length(osiris_o3_bk.tot_col)
        osiris_o3_bk.num_dens2(i,:)=interp1(14.5:51.5,osiris_o3_bk.num_dens(i,:),alt_bk(20:39));
%         osiris_o3_bk.num_dens2(i,:)=apriori(i,20:39)+avk(i,20:39).*(osiris_o3_bk.num_dens2(i,:)-apriori(i,20:39));
    end
    osiris_o3_bk.num_dens=[];
    osiris_o3_bk.Properties.VariableNames{'num_dens2'} = 'num_dens';

    %%%%
    [alt_bk,~,apriori,avk,part_prof,~]=read_bruker_prof_avk(1,ace_fts_o3_bk.mjd2k,'bruker');
    ace_fts_o3_bk.num_dens2=NaN(length(ace_fts_o3_bk.tot_col),20);
    for i=1:length(ace_fts_o3_bk.tot_col)
        ace_fts_o3_bk.num_dens2(i,:)=interp1(14.5:51.5,ace_fts_o3_bk.num_dens(i,:),alt_bk(20:39));
%         ace_fts_o3_bk.num_dens2(i,:)=apriori(i,20:39)+avk(i,20:39).*(ace_fts_o3_bk.num_dens2(i,:)-apriori(i,20:39));
    end
    ace_fts_o3_bk.num_dens=[];
    ace_fts_o3_bk.Properties.VariableNames{'num_dens2'} = 'num_dens';
    
    %%%%
    [alt_bk,~,apriori,avk,part_prof,~]=read_bruker_prof_avk(1,ace_mae_o3_bk.mjd2k,'bruker');
    ace_mae_o3_bk.num_dens2=NaN(length(ace_mae_o3_bk.tot_col),20);
    for i=1:length(ace_mae_o3_bk.tot_col)
        ace_mae_o3_bk.num_dens2(i,:)=interp1(14.5:51.5,ace_mae_o3_bk.num_dens(i,:),alt_bk(20:39));
%         ace_mae_o3_bk.num_dens2(i,:)=apriori(i,20:39)+avk(i,20:39).*(ace_mae_o3_bk.num_dens2(i,:)-apriori(i,20:39));
    end
    ace_mae_o3_bk.num_dens=[];
    ace_mae_o3_bk.Properties.VariableNames{'num_dens2'} = 'num_dens';

    figure()
    
    if save_as_pdf
        fig_ax = tight_subplot(3,3,[0.07,0.0],[0.08,0.02],[0.13,0.07]);
    else
        fig_ax = tight_subplot(3,3,[0.07,0.0],[0.10,0.04],[0.13,0.07]);
    end

    ii=1;
    diff_prof_time(fig_ax(1:3),osiris_o3_bk,bruker_o3,label_osiris,label_bruker,...
                       'left', 0, alt_bk(20:39))

    ii=2;
    diff_prof_time(fig_ax(4:6),ace_fts_o3_bk,bruker_o3,label_fts,label_bruker,...
                       'left', 0, alt_bk(20:39))

    ii=3;
    diff_prof_time(fig_ax(7:9),ace_mae_o3_bk,bruker_o3,label_mae,label_bruker,...
                       'bottomleft', 0, alt_bk(20:39))

    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 850]);

    save_pdf(save_figs,'diff_prof_bk_o3')


%%%     select MAESTRO profiles with large errors
%%%
% %     occ=[51475,40549,3123,49244,67799,65312,13645,56816,35390,51506];
% %     [~,ind,~]=intersect(ace_mae_o3.occultation,occ);
% %     ace_mae_o3=ace_mae_o3(ind,:);
% %     
% %     mae=[ace_mae_o3{:,1:3}, ace_mae_o3{:,7}];
% %     fts=[ace_fts_o3{:,1:3}, ace_fts_o3{:,7}];
% %     [~,ind,~]=intersect(fts,mae,'rows');
% %     ace_fts_o3=ace_fts_o3(ind,:);
% % 
% %     alt=14.5:51.5;
% %     for i=1:10
% %         figure()
% %         subplot(2,2,[1 3])
% %         plot(ace_fts_o3.num_dens(i,:),alt,'b-','linewidth',1.3), hold on
% %         plot(ace_mae_o3.num_dens(i,:),alt,'r-','linewidth',1.3), hold on
% %         
% %         legend('ACE-FTS','MAESTRO')
% %         xlabel('O_3 (molec/cm^3)')
% %         ylabel('Altitude (km)')
% %         ylim([10,55])
% %         xlim([0,8]*1e12)
% %         grid on
% %         ampm='sr';
% %         if ace_mae_o3.ampm(i)==1, ampm='ss'; end
% %         title(['Profile: ' ampm num2str(ace_mae_o3.occultation(i))])
% % 
% %         subplot(222)
% %         plot(ace_fts_o3.num_dens(i,:)-ace_mae_o3.num_dens(i,:),alt,'k-','linewidth',1.3)
% %         xlabel('\Delta_{abs}')
% %         ylim([10,55])
% %         xlim([-0.6,2.2]*1e12)
% %         grid on
% %         
% %         subplot(224)
% %         plot(200*(ace_fts_o3.num_dens(i,:)-ace_mae_o3.num_dens(i,:))...
% %              ./(ace_fts_o3.num_dens(i,:)+ace_mae_o3.num_dens(i,:)),alt,'k-','linewidth',1.3)
% %         xlabel('\Delta_{rel} (%)')
% %         ylim([10,55])
% %         xlim([-70,70])
% %         grid on
% %         
% %         saveas(gcf,[ampm num2str(ace_mae_o3.occultation(i)) '.png'])
% %         
% %     end

  
end

%% no2 gb vs sat total column correlation plots
if no2_corr
    
    figure()
    
    fig_ax = tight_subplot(2,4,[0.09,0.042],[0.11,0.075],[0.095,0.05]);

    ii=1;
    
    % OSIRIS vs GB
    axes(fig_ax(ii))
    corr_time(gbs_no2,osiris_no2,label_gbs,label_osiris,'left',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(gbs_no2uv,osiris_no2,label_gbs_uv,label_osiris,'',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(saoz_no2,osiris_no2,label_saoz,label_osiris,'',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_time(bruker_no2,osiris_no2_bk,label_bruker,label_osiris,'',false);
    ii=ii+1;

    % ACE-FTS vs GB
    axes(fig_ax(ii))
    corr_twilight(gbs_no2,ace_fts_no2,label_gbs,label_fts,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_twilight(gbs_no2uv,ace_fts_no2,label_gbs_uv,label_fts,'bottom',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_twilight(saoz_no2,ace_fts_no2,label_saoz,label_fts,'bottom',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(bruker_no2,ace_fts_no2_bk,label_bruker,label_fts,'bottom',false);
    ii=ii+1;
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 520]);
    
    save_pdf(save_figs,'corr_no2')

end


if no2_diff
    
    figure()
    
    fig_ax = tight_subplot(2,2,[0.16,0.05],[0.11,0.12],[0.1,0.05]);

    ii=1;
    
    % OSIRIS vs GB
    axes(fig_ax(ii))
    diff_time(gbs_no2,osiris_no2,label_gbs,label_osiris,'left',false,1);

    axes(fig_ax(ii))
    diff_time(gbs_no2uv,osiris_no2,label_gbs_uv,label_osiris,'',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    diff_time(saoz_no2,osiris_no2,label_saoz,label_osiris,'',false,1);
    
    axes(fig_ax(ii))
    diff_time(bruker_no2,osiris_no2_bk,label_bruker,label_osiris,'',false,2);
    ii=ii+1;

    % ACE-FTS vs GB
    axes(fig_ax(ii))
    diff_twilight(gbs_no2,ace_fts_no2,label_gbs,label_fts,'bottomleft',false,1);

    axes(fig_ax(ii))
    diff_twilight(gbs_no2uv,ace_fts_no2,label_gbs_uv,label_fts,'bottom',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    diff_twilight(saoz_no2,ace_fts_no2,label_saoz,label_fts,'bottom',false,1);

    axes(fig_ax(ii))
    diff_time(bruker_no2,ace_fts_no2_bk,label_bruker,label_fts,'',false,2);
    ii=ii+1;
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 550]);
    
    save_pdf(save_figs,'diff_no2')
    
end


if smooth_diffs
    
    load('/home/kristof/work/satellite_validation/diff_tables.mat')
    
    % o3 changes
    header=o3_rel_diff_nosmooth.Properties.VariableNames;

    % get labels, add minus sign
    for i=1:length(header)
        header{i}=[header{i}(1:2) '-' header{i}(3:4)];
    end
    
    figure()
    subplot(211)
    errorbar(1:length(header),o3_rel_diff_nosmooth{1,:},o3_rel_diff_nosmooth{2,:},'k.'), hold on
    errorbar(1:length(header),o3_rel_diff_smooth{1,:},o3_rel_diff_smooth{2,:},'r.')
    
    set(gca,'Xtick',1:length(header),'XTickLabel',header,'XTickLabelRotation',90,'YGrid','on')
    ylabel('\Delta_{rel} (%)')
    
    set(gca, 'FontSize', label_size+2)
    
    % no2 changes
    header=no2_rel_diff_nosmooth.Properties.VariableNames;

    % get labels, add minus sign
    for i=1:length(header)
        header{i}=[header{i}(1:2) '-' header{i}(3:4)];
    end
    
    subplot(212)
    errorbar(1:length(header),no2_rel_diff_nosmooth{1,:},no2_rel_diff_nosmooth{2,:},'k.'), hold on
    errorbar(1:length(header),no2_rel_diff_smooth{1,:},no2_rel_diff_smooth{2,:},'r.')
    
    set(gca,'Xtick',1:length(header),'XTickLabel',header,'XTickLabelRotation',90,'YGrid','on')
    ylabel('\Delta_{rel} (%)')
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gca, 'FontSize', label_size+2)
    set(gcf, 'Position', [100, 100, 520, 500]);

% %     % cloud filter changes for ozone
% %     header=o3_rel_diff_noCF.Properties.VariableNames;
% % 
% %     % get labels, add minus sign
% %     for i=1:length(header)
% %         header{i}=[header{i}(1:2) '-' header{i}(3:4)];
% %     end
% %     
% %     figure()
% %     errorbar(1:length(header),o3_rel_diff_noCF{1,:},o3_rel_diff_noCF{2,:},'k.'), hold on
% %     errorbar(1:length(header),o3_rel_diff_CF{1,:},o3_rel_diff_CF{2,:},'r.')
    
    
end

if tc_stuff
    
    figure()
    axis off
    
    % indices of tc fields, col1: row numbers for sat and gb; row1/2: sat-gb-x triplet
    % indices for sat/gb
    %%% OSIRIS ozone
    tmp=[[6,1,2,3,4,12,13];...
         [1,14,17,20,23,24,25]];
    mean_tc(1,tmp,'OS, GV (O3)',1);

    tmp=[[6,1,5,6,7,15,16];...
         [2,11,17,20,23,24,25]];
    mean_tc(1,tmp,'OS, SA (O3)',2);

    tmp=[[6,2,5,8,9,18,19];...
         [3,11,14,20,23,24,25]];
    mean_tc(1,tmp,'OS, BK (O3)',3);

    tmp=[[6,3,6,8,10,21,22];...
         [4,11,14,17,23,24,25]];
    mean_tc(1,tmp,'OS, PA (O3)',4);

    tmp=[[6,4,7,9,10];...
         [5,11,14,17,20]];
    mean_tc(1,tmp,'OS, BW (O3)',5);

    %%% ACE-FTS ozone
    tmp=[[7,1,2,3,11,13];...
         [1,15,18,21,24,26]];
    mean_tc(1,tmp,'AF, GV (O3)',6);
    
    tmp=[[7,1,5,6,14,16];...
         [2,12,18,21,24,26]];
    mean_tc(1,tmp,'AF, SA (O3)',7);
    
    tmp=[[7,2,5,8,17,19];...
         [3,12,15,21,24,26]];
    mean_tc(1,tmp,'AF, BK (O3)',8);

    tmp=[[7,3,6,8,20,22];...
         [4,12,15,18,24,26]];
    mean_tc(1,tmp,'AF, PA (O3)',9);

    %%% ACE-MAESTRO ozone
    tmp=[[8,1,2,3,11,12];...
         [1,16,19,22,25,26]];
    mean_tc(1,tmp,'AM, GV (O3)',10);
    
    tmp=[[8,1,5,6,14,15];...
         [2,13,19,22,25,26]];
    mean_tc(1,tmp,'AM, SA (O3)',11);
    
    tmp=[[8,2,5,8,17,18];...
         [3,13,16,22,25,26]];
    mean_tc(1,tmp,'AM, BK (O3)',12);
    
    tmp=[[8,3,6,8,20,21];...
         [4,13,16,19,25,26]];
    mean_tc(1,tmp,'AM, PA (O3)',13);

    %%% OSIRIS NO2
    tmp=[[5,1,2,3];...
         [1,9,11,13]];
    mean_tc(2,tmp,'OS, GV (NO2)',14);

    tmp=[[5,1,4,5];...
         [2,7,11,13]];
    mean_tc(2,tmp,'OS, GU (NO2)',15);
    
    tmp=[[5,2,4,6];...
         [3,7,9,13]];
    mean_tc(2,tmp,'OS, SA (NO2)',16);

    tmp=[[5,3,5,6];...
         [4,7,9,11]];
    mean_tc(2,tmp,'OS, BK (NO2)',17);

    %%% ACE-FTS NO2
    tmp=[[6,1,2,3];...
         [1,10,12,14]];
    mean_tc(2,tmp,'AF, GV (NO2)',18);
    
    tmp=[[6,1,4,5];...
         [2,8,12,14]];
    mean_tc(2,tmp,'AF, GU (NO2)',19);
    
    tmp=[[6,2,4,6];...
         [3,8,10,14]];
    mean_tc(2,tmp,'AF, SA (NO2)',20);

    tmp=[[6,3,5,6];...
         [3,8,10,12]];
    mean_tc(2,tmp,'AF, BK (NO2)',21);
    
    %%% Sat sat O3
    tmp=[[6,12,15,18,21,26];...
         [7,11,14,17,20,25]];
    mean_tc(1,tmp,'OS, AF (O3)',22);

    tmp=[[6,13,16,19,22,26];...
         [8,11,14,17,20,24]];
    mean_tc(1,tmp,'OS, AM (O3)',23);

    tmp=[[7,13,16,19,22,25];...
         [8,12,15,18,21,24]];
    mean_tc(1,tmp,'AF, AM (O3)',24);
    
    title('TCA with given pair + all possible 3rd instruments')
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 720, 800]);
    
    save_pdf(save_figs,'TCA_rmse_rt')    
    
end


%% ozone gb total column plots
if o3_gb
    
    %%% correlation plots
    figure()
    
%     fig_ax = tight_subplot(2,3,[0.14,0.035],[0.11,0.09],[0.12,0.08]); no PA
    fig_ax = tight_subplot(3,4,[0.04,0.06],[0.04,0.03],[0.07,0.03]);

    axes(fig_ax(8))
    delete(gca)
    axes(fig_ax(12))
    delete(gca)
    
    ii=1;

    % all minus brewer
    axes(fig_ax(ii))
    corr_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(brewer_o3_ds,paris_o3,label_brewer,label_paris,'bottomleft',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'bottomleft',false);
    ii=ii+1;
    
    
    % all minus bruker
    axes(fig_ax(ii))
    corr_time(bruker_o3,gbs_o3,label_bruker,label_gbs,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(bruker_o3,saoz_o3,label_bruker,label_saoz,'bottomleft',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_time(bruker_o3,paris_o3,label_bruker,label_paris,'bottomleft',bruker_part_col);
    ii=ii+1;

    % all minus PA
    axes(fig_ax(ii+1))
    corr_time(paris_o3,gbs_o3,label_paris,label_gbs,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii+1))
    corr_time(paris_o3,saoz_o3,label_paris,label_saoz,'bottomleft',false);
    ii=ii+1;
    
    % doas
    axes(fig_ax(ii+1))
    corr_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottomleft',false);

    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
%     set(gcf, 'Position', [100, 100, 1000, 550]);no PA
    set(gcf, 'Position', [100, 100, 1000, 800]);
    
    save_pdf(save_figs,'gb_corr_o3')    

    
    %%% difference plots
    figure()
    
%     fig_ax = tight_subplot(3,1,[0.11,0.05],[0.04,0.05],[0.16,0.08]); no PA
    
    if save_as_pdf, 
        fig_ax = tight_subplot(3,2,[0.11,0.05],[0.04,0.05],[0.1,0.05]);    
    else
        fig_ax=tight_subplot(3,2,[0.11,0.05],[0.06,0.07],[0.1,0.05]);
    end
    
    ii=1;
    
    % all minus brewer
    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,gbs_o3,label_brewer,label_gbs,'',false,1);

    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,saoz_o3,label_brewer,label_saoz,'left',false,2);
    ii=ii+1;

    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,paris_o3,label_brewer,label_paris,'',false,1);

    axes(fig_ax(ii))
    diff_time(brewer_o3_ds,bruker_o3,label_brewer,label_bruker,'',false,2);
    ii=ii+1;
    
    % all minus bruker
    axes(fig_ax(ii))
    diff_time(bruker_o3,gbs_o3,label_bruker,label_gbs,'',false,1);

    axes(fig_ax(ii))
    diff_time(bruker_o3,saoz_o3,label_bruker,label_saoz,'left',false,2);
    ii=ii+1;
    
    axes(fig_ax(ii))
    diff_time(bruker_o3,paris_o3,label_bruker,label_paris,'',bruker_part_col,11);
    ii=ii+1;
    
    % all minus PARIS
    axes(fig_ax(ii))
    diff_time(paris_o3,gbs_o3,label_paris,label_gbs,'',false,1);

    axes(fig_ax(ii))
    diff_time(paris_o3,saoz_o3,label_paris,label_saoz,'bottomleft',false,2);
    ii=ii+1;
    
    % doas
    axes(fig_ax(ii))
    diff_twilight(saoz_o3,gbs_o3,label_saoz,label_gbs,'bottom',false,11);

    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
%     set(gcf, 'Position', [100, 100, 520, 850]); no PA
    set(gcf, 'Position', [100, 100, 1000, 850]);
    
    save_pdf(save_figs,'gb_diff_o3')


end

%% no2 gb total column plots
if no2_gb
    
    %%% correlation plots
    figure()
    
    fig_ax = tight_subplot(2,3,[0.14,0.035],[0.11,0.09],[0.12,0.08]);

    ii=1;
    
    % SAOZ vs GBS
    axes(fig_ax(ii))
    corr_twilight(saoz_no2,gbs_no2,label_saoz,label_gbs,'bottomleft',false);
    ii=ii+1;

    % bruker vs GBS
    axes(fig_ax(ii))
    corr_time(gbs_no2,bruker_no2,label_gbs,label_bruker,'bottomleft',false);
    ii=ii+1;

    % other
    axes(fig_ax(ii))
    corr_twilight(gbs_no2,gbs_no2uv,label_gbs,label_gbs_uv,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_twilight(saoz_no2,gbs_no2uv,label_saoz,label_gbs_uv,'bottomleft',false);
    ii=ii+1;
    
    axes(fig_ax(ii))
    corr_time(gbs_no2uv,bruker_no2,label_gbs_uv,label_bruker,'bottomleft',false);
    ii=ii+1;

    axes(fig_ax(ii))
    corr_time(saoz_no2,bruker_no2,label_saoz,label_bruker,'bottomleft',false);
    ii=ii+1;
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 1000, 550]);
    
    save_pdf(save_figs,'gb_corr_no2')    
   
    
    %%% difference plots
    figure()
    
    if save_as_pdf, 
        fig_ax = tight_subplot(3,1,[0.11,0.05],[0.04,0.05],[0.16,0.08]);
    else
        fig_ax=tight_subplot(3,1,[0.11,0.05],[0.06,0.07],[0.16,0.08]);
    end

    ii=1;
    
    % SAOZ vs GBS
    axes(fig_ax(ii))
    diff_twilight(saoz_no2,gbs_no2,label_saoz,label_gbs,'',false,1);

    axes(fig_ax(ii))
    diff_twilight(saoz_no2,gbs_no2uv,label_saoz,label_gbs_uv,'left',false,2);
    ii=ii+1;
    ylim([-30,30])
    
    % bruker vs GBS
    axes(fig_ax(ii))
    diff_time(gbs_no2,bruker_no2,label_gbs,label_bruker,'',false,1);

    axes(fig_ax(ii))
    diff_time(gbs_no2uv,bruker_no2,label_gbs_uv,label_bruker,'left',false,2);
    ii=ii+1;
    ylim([-30,30])

    % other
    axes(fig_ax(ii))
    diff_twilight(gbs_no2,gbs_no2uv,label_gbs,label_gbs_uv,'',false,1);

    axes(fig_ax(ii))
    diff_time(saoz_no2,bruker_no2,label_saoz,label_bruker,'bottomleft',false,2);
    ii=ii+1;
    ylim([-30,30])
    
    
    set(findall(gcf,'-property','FontName'),'FontName','Arial')
    set(gcf, 'Position', [100, 100, 520, 850]);
    
    save_pdf(save_figs,'gb_diff_no2')

    
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% local functions to avoid copy-pasting code
function corr_twilight(table1, table2, instr1, instr2, loc_str, partcol)

    global alphabet text_size label_size ii

    [ x1, x2, e1, e2, times, DMP ] =...
        find_coincidences_twilight( table1, table2, partcol );

%     T2=DMP.T2(:,13);
%     spv1=DMP.spv1(:,13);
%     lat1=DMP.lat1(:,13);
%     lat2=DMP.lat2;
%     
%     ind=[];
% %     ind=find(spv1>1.2e-4 & spv1<1.6e-4);
% %     ind=find(abs(lat1-lat2)>2);
%     x1(ind)=[];
%     x2(ind)=[];
%     e1(ind)=[];
%     e2(ind)=[];
    
    title_str=[alphabet(ii) ')'];

    corr_plot([x1,e1],[x2,e2],instr1, instr2, loc_str, partcol,title_str,text_size,label_size);
    
end

function corr_time(table1, table2, instr1, instr2, loc_str, partcol)

    global alphabet text_size label_size ii

    dt=12;
    
    [ x1, x2, e1, e2, times,~, DMP ] =...
        find_coincidences_time( table1, table2, dt, partcol );

%     T2=DMP.T2(:,13);
%     spv1=DMP.spv1(:,13);
%     lat1=DMP.lat1(:,13);
%     lat2=DMP.lat2;
% 
%     ind=[];
% %     ind=find(spv1>1.2e-4 & spv1<1.6e-4);
% %     ind=find(abs(lat1-lat2)>2);
%     x1(ind)=[];
%     x2(ind)=[];
%     e1(ind)=[];
%     e2(ind)=[];

    title_str=[alphabet(ii) ')'];
    
    corr_plot([x1,e1],[x2,e2],instr1, instr2, loc_str, partcol,title_str,text_size,label_size);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function diff_twilight(table1, table2, instr1, instr2, loc_str, partcol, double_plot)

    global alphabet text_size label_size ii dmp_filter
    
    [ x1, x2, e1, e2, times, DMP ] =...
        find_coincidences_twilight( table2, table1, partcol );
    
    if dmp_filter==1
        good_ind=dmp_criteria(DMP, instr1, instr2);
        x1=x1(good_ind);
        x2=x2(good_ind);
        times=times(good_ind);
    elseif dmp_filter==2
        good_ind=lat_filter(DMP, instr1, instr2);
        x1=x1(good_ind);
        x2=x2(good_ind);
        times=times(good_ind);
    end
    
    sza1=DMP.sza1; % satellite SZA at tangent point
    sza2=DMP.sza2; % ground-based SZA
    
    title_str=[alphabet(ii) ')'];
    
%     diff_plot_time('abs',x1,x2,times,'timeseries',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
    diff_plot_time('abs',x1,x2,times,'seasonal',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
%     diff_plot_time('abs',x1,x2,sza2,'SZA',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)    
%     diff_plot_time('rel',x1,x2,times,'trend',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)

end

function diff_time(table1, table2, instr1, instr2, loc_str, partcol, double_plot)

    global alphabet text_size label_size ii dmp_filter
    
    dt=12;
    
    [ x1, x2, e1, e2, times,~, DMP ] =...
        find_coincidences_time( table2, table1, dt, partcol );
  
    if dmp_filter==1
        good_ind=dmp_criteria(DMP, instr1, instr2);
        x1=x1(good_ind);
        x2=x2(good_ind);
        times=times(good_ind);
    elseif dmp_filter==2
        good_ind=lat_filter(DMP, instr1, instr2);
        x1=x1(good_ind);
        x2=x2(good_ind);
        times=times(good_ind);
    end
    
    sza1=DMP.sza1; % satellite SZA at tangent point
    sza2=DMP.sza2; % ground-based SZA
    
    title_str=[alphabet(ii) ')'];
    
%     diff_plot_time('abs',x1,x2,times,'timeseries',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
    diff_plot_time('abs',x1,x2,times,'seasonal',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
%     diff_plot_time('abs',x1,x2,sza2,'SZA',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
%     diff_plot_time('rel',x1,x2,times,'trend',instr2,instr1,loc_str,title_str,text_size,label_size,double_plot)
     
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function diff_prof_twilight(axes, table1, table2, instr1, instr2, loc_str, ace, double_plot)

    global alphabet text_size label_size ii
    
    [ ~, ~, ~, ~, ~, ~, prof1, prof2 ] =...
        find_coincidences_twilight( table1, table2, false, ace );
    
    title_str=[alphabet(ii) ')'];
    
    disp([instr1 '-' instr2])
    
    diff_plot_prof( axes, prof1,prof2, instr1,instr2,...
                         loc_str,title_str,text_size,label_size,double_plot);
    

end


function diff_prof_time(axes, table1, table2, instr1, instr2, loc_str, double_plot,alt)

    global alphabet text_size label_size ii
    
    [ ~, ~, ~, ~, ~, ~, ~, prof1, prof2 ] =...
        find_coincidences_time( table1, table2, 12, false );
    
    title_str=[alphabet(ii) ')'];
    
    disp([instr1 '-' instr2])
    
    diff_plot_prof( axes, prof1,prof2, instr1,instr2,...
                         loc_str,title_str,text_size,label_size,double_plot,alt);
    

end

function save_pdf(save_figs, fname)
    global save_as_pdf

    if save_figs
        
        h=gcf;

        set(h,'Units','Inches');

        pos = get(h,'Position');

        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

        if save_as_pdf
            % pdf images
            f_out=['/home/kristof/work/documents/paper_sat-val/figures/' fname '.pdf'];
            print(h,f_out,'-dpdf','-r0')
        
        else
            % jpg images
            f_out=['/home/kristof/work/documents/paper_sat-val/figures/' fname '.jpg'];
            print(h,f_out,'-djpeg','-r400','-opengl')
%             % png images
%             f_out=['/home/kristof/work/documents/paper_sat-val/figures/' fname '.png'];
%             print(h,f_out,'-dpng','-r400','-opengl')
            
        end
        
        close(gcf)
        
    else
        return
    end

end

function [rmse,rt,rmse_err,rt_err,N_sum]=mean_tc(tg,ind,pair_str,num)

    global text_size
    
    data_file='/home/kristof/work/satellite_validation/all_data_tc.mat';
    load(data_file)

    % mean RMSE and R^t and total number of coincidences
    rmse=zeros(2,1);
    rt=zeros(2,1);
    N_sum=0;
    N_coinc=size(ind,2)-1;

    % standard error
    rmse_err=zeros(2,1);
    rt_err=zeros(2,1);

    for i=1:2
        if tg==1
            rmse(i)=nanmean(tc_o3_rmse{ind(i,1),ind(i,2:end)});
            rt(i)=nanmean(tc_o3_rt{ind(i,1),ind(i,2:end)});
            
            N_sum=nansum(tc_o3_N{ind(i,1),ind(i,2:end)});

            rmse_err(i)=std_err(tc_o3_rmse{ind(i,1),ind(i,2:end)});
            rt_err(i)=std_err(tc_o3_rt{ind(i,1),ind(i,2:end)});

        else
            rmse(i)=nanmean(tc_no2_rmse{ind(i,1),ind(i,2:end)});
            rt(i)=nanmean(tc_no2_rt{ind(i,1),ind(i,2:end)});
            
            N_sum=nansum(tc_no2_N{ind(i,1),ind(i,2:end)});

            rmse_err(i)=std_err(tc_no2_rmse{ind(i,1),ind(i,2:end)});
            rt_err(i)=std_err(tc_no2_rt{ind(i,1),ind(i,2:end)});

        end
    end
    
    text(0,1-0.045*num,sprintf([pair_str...
        ':   sum(RMSE) = %.3g +- %.2g;   R_t= %.2f+-%.2f, %.2f+-%.2f;  N=%.0f (%.0f)\n'],...
        [sum(rmse),sqrt(sum(rmse_err.^2)),rt(1),rt_err(1),rt(2),rt_err(2),N_sum,N_coinc])...
        ,'Units','normalized','fontsize',text_size)
    
    
end

function good_ind=dmp_criteria(DMP, instr1, instr2)

    global alt_dmp_in
    
    % initialize index of good values as index of all values
    good_ind=1:size(DMP.T1,1);
    
    % coincidence criteria
    
    t_lim=10; % max temperature difference
    spv_lim=[1.2,1.6]*1e-4; % vortex boundary sPV: either both in or both out
    
    % altitudes where coincidences are checked; different for ACE and
    % OSIRIS -- return if not comparig to satellite
    if strcmp(instr2,'AF') || strcmp(instr2,'AM')
        alt_list=[14,18,20,22];
%         alt_list=[24,26,30];

%         if strcmp(instr1,'BK')
%             alt_list=[14,18,20,22,24,26,30,36,46]; 
%         else
%             alt_list=[14,18,20,22]; 
%         end
    elseif strcmp(instr2,'OS')
        alt_list=[25];
    else
        return
    end
    
    % return if comparing to brewer -- no DMPs
    if strcmp(instr1,'BW'), return, end
    
    % loop over altitudes
    for i=alt_list
        
        ind=find(alt_dmp_in==i);
        
        good_ind_tmp=find( ( (DMP.spv1(:,ind)<spv_lim(1) & DMP.spv2(:,ind)<spv_lim(1)) | ...
                             (DMP.spv1(:,ind)>spv_lim(2) & DMP.spv2(:,ind)>spv_lim(2))      ) & ...
                           abs(DMP.T1(:,ind)-DMP.T2(:,ind))<t_lim );
        
        good_ind=intersect(good_ind,good_ind_tmp);
        
    end
    
end

function good_ind=lat_filter(DMP, instr1, instr2)

    global alt_dmp_in
    
    % initialize index of good values as index of all values
    good_ind=1:size(DMP.T1,1);
    
    if strcmp(instr1,'BW'), return, end

    if strcmp(instr2,'AF') || strcmp(instr2,'AM') 
        alt_list=[30]; 
    elseif strcmp(instr2,'OS')
        alt_list=[25];
    else
        return
    end
    
    for i=alt_list
        
        ind=find(alt_dmp_in==i);
        
%         good_ind_tmp=find(abs(DMP.lat1(:,ind)-DMP.lat2(:,ind))<2);
        good_ind_tmp=find(abs(DMP.lat1-DMP.lat2(:,ind))<2);
        good_ind=intersect(good_ind,good_ind_tmp);
    end
end

