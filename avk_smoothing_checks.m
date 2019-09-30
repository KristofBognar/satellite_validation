

ozone=0;

if ozone
    
    load('/home/kristof/work/satellite_validation/all_data_nosmooth.mat')
    AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';

    % load ozonesonde data for DOAS LUTs
    load('/home/kristof/work/ozonesonde/Eureka/sonde_for_VCD.mat')
    % time when we have the sonde measurement
    sonde_time = sonde(:,1) +(sonde(:,2)-1+sonde(:,3)/24)./daysinyear(sonde(:,1));
    sonde_ozone = sonde(:,4)*2.69e16;% the sonde ozone VCD (in molec/cm2)

    % AF measurement time
    measurement_time = ace_fts_o3.year + ace_fts_o3.fractional_time./daysinyear(ace_fts_o3.year);
    sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

    % AF smoothing profiles and AVKs
    [lut_prof_af,lut_avk_af]=read_DOAS_prof_avk(1,[ace_fts_o3.year,ace_fts_o3.day,sonde_in],...
                                          AVK_LUT_dir );

    % AM measurement time
    measurement_time = ace_mae_o3.year + ace_mae_o3.fractional_time./daysinyear(ace_mae_o3.year);
    sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

    % AM smoothing profiles and AVKs
    [lut_prof_am,lut_avk_am]=read_DOAS_prof_avk(1,[ace_mae_o3.year,ace_mae_o3.day,sonde_in],...
                                          AVK_LUT_dir );


    % OSIIRIS (use spring/fall only)                      
    day_range=[106.25,239.25];

    ind=find(osiris_o3.fractional_time>day_range(1)-1 & osiris_o3.fractional_time<day_range(2));
    osiris_o3(ind,:)=[];

    % OS measurement time
    measurement_time = osiris_o3.year + osiris_o3.fractional_time./daysinyear(osiris_o3.year);
    sonde_in = interp1(sonde_time, sonde_ozone, measurement_time);

    % OS smoothing profiles and AVKs
    [lut_prof_os,lut_avk_os]=read_DOAS_prof_avk(1,[osiris_o3.year,...
                                          floor(osiris_o3.fractional_time)+1,...
                                          sonde_in], AVK_LUT_dir );

    alt=0.5:59.5;

    figure()
    subplot(121)
    plot(mean(lut_prof_af(:,4:end)),alt,'b-'), hold on
    plot(mean(lut_prof_os(:,4:end)),alt,'r-'), hold on
    plot(mean(lut_prof_am(:,4:end)),alt,'g--'), hold on

    plot(mean(ace_fts_o3.num_dens),14.5:51.5,'bx-')
    plot(mean(osiris_o3.num_dens),14.5:51.5,'rx-')
    plot(mean(ace_mae_o3.num_dens),14.5:51.5,'gx-')

    legend('AF mean a priori','OS mean a priori','AM mean a priori',...
           'AF mean profile','OS mean profile','AM mean profile')
    xlabel('O_3 num dens (molec/cm^3)')
    ylabel('Altitude (km)')
    ylim([10,52])
    grid on

    subplot(122)
    plot(mean(lut_avk_af(:,4:end)),alt,'b-'), hold on
    plot(mean(lut_avk_os(:,4:end)),alt,'r-'), hold on
    plot(mean(lut_avk_am(:,4:end)),alt,'g--'), hold on

    legend('AF','OS','AM')
    xlabel('Column AVK')
    ylim([10,52])
    grid on



    mean_ap_af=mean(lut_prof_af(:,18:end-8));
    mean_ap_am=mean(lut_prof_am(:,18:end-8));
    mean_ap_os=mean(lut_prof_os(:,18:end-8));
    mean_avk_af=mean(lut_avk_af(:,18:end-8));
    mean_avk_am=mean(lut_avk_am(:,18:end-8));
    mean_avk_os=mean(lut_avk_os(:,18:end-8));

    mean_af=mean(ace_fts_o3.num_dens);
    mean_am=mean(ace_mae_o3.num_dens);
    mean_os=mean(osiris_o3.num_dens);

    figure()
    new_af=mean_ap_af+mean_avk_af.*(mean_af-mean_ap_af);
    plot(mean(lut_prof_af(:,4:end)),alt,'b-'), hold on
    plot(mean(ace_fts_o3.num_dens),14.5:51.5,'bx-')
    plot(new_af,14.5:51.5,'kx-')
    title('AF mean smoothed')
    legend('mean a priori','mean profile','smoothed')

    figure()
    new_am=mean_ap_am+mean_avk_am.*(mean_am-mean_ap_am);
    plot(mean(lut_prof_am(:,4:end)),alt,'b-'), hold on
    plot(mean(ace_mae_o3.num_dens),14.5:51.5,'bx-')
    plot(new_am,14.5:51.5,'kx-')
    title('AM mean smoothed')
    legend('mean a priori','mean profile','smoothed')

    figure()
    new_os=mean_ap_os+mean_avk_os.*(mean_os-mean_ap_os);
    plot(mean(lut_prof_os(:,4:end)),alt,'b-'), hold on
    plot(mean(osiris_o3.num_dens),14.5:51.5,'bx-')
    plot(new_os,14.5:51.5,'kx-')
    title('OS mean smoothed')
    legend('mean a priori','mean profile','smoothed')

else % no2
    
%     load('all_data_nosmooth.mat')
%     AVK_LUT_dir='/home/kristof/work/NDACC/guidelines/2012/';
% 
%     %%% AF AVK smoothing
%     [lut_prof_af,lut_avk_af]=read_DOAS_prof_avk(2,[ace_fts_no2.year,ace_fts_no2.fractional_time],...
%                                           AVK_LUT_dir );
% 
%     % same for UV
%     [lut_prof_afuv,lut_avk_afuv]=read_DOAS_prof_avk(3,[ace_fts_no2.year,ace_fts_no2.fractional_time],...
%                                           AVK_LUT_dir );
% 
% 
%     % OSIIRIS (use spring/fall only)                      
%     day_range=[106.25,239.25];
% 
%     ind=find(osiris_no2.fractional_time>day_range(1)-1 & osiris_no2.fractional_time<day_range(2));
%     osiris_no2(ind,:)=[];
% 
%     %%% OS AVK smoothing
%     [lut_prof_os,lut_avk_os]=read_DOAS_prof_avk(2,[osiris_no2.year,osiris_no2.fractional_time],...
%                                           AVK_LUT_dir );
% 
%     % same for UV
%     [lut_prof_osuv,lut_avk_osuv]=read_DOAS_prof_avk(3,[osiris_no2.year,osiris_no2.fractional_time],...
%                                           AVK_LUT_dir );
% 
%     alt=0.5:59.5;

    figure()
    subplot(121)
    % plot vis a priori only, LUT profile is not wavelength dependent
    plot(mean(lut_prof_af(:,3:end)),alt,'b-'), hold on
    plot(mean(lut_prof_os(:,3:end)),alt,'r-'), hold on

    plot(mean(ace_fts_no2.num_dens),12.5:39.5,'bx-')
    plot(mean(osiris_no2.num_dens),12.5:31.5,'rx-')

    legend('AF mean a priori','OS mean a priori',...
           'AF mean profile','OS mean profile')
    xlabel('NO_2 num dens (molec/cm^3)')
    ylabel('Altitude (km)')
    ylim([10,52])
    grid on

    subplot(122)
    % LUT AVK is wavelength dependent
    plot(mean(lut_avk_af(:,3:end)),alt,'b-'), hold on
    plot(mean(lut_avk_afuv(:,3:end)),alt,'b--'), hold on
    plot(mean(lut_avk_os(:,3:end)),alt,'r-'), hold on
    plot(mean(lut_avk_osuv(:,3:end)),alt,'r--'), hold on

    legend('AF','AFuv','OS','OSuv')
    xlabel('Column AVK')
    ylim([10,52])
    grid on



%     mean_ap_af=mean(lut_prof_af(:,13:end-20));
%     mean_ap_os=mean(lut_prof_os(:,13:end-28));
%     mean_avk_af=mean(lut_avk_af(:,13:end-20));
%     mean_avk_os=mean(lut_avk_os(:,13:end-28));
% 
%     mean_af=mean(ace_fts_no2.num_dens);
%     mean_os=mean(osiris_no2.num_dens);
% 
%     figure()
%     new_af=mean_ap_af+mean_avk_af.*(mean_af-mean_ap_af);
%     plot(mean(lut_prof_af(:,3:end)),alt,'b-'), hold on
%     plot(mean(ace_fts_no2.num_dens),12.5:39.5,'bx-')
%     plot(new_af,12.5:39.5,'kx-')
%     title('AF mean smoothed')
%     legend('mean a priori','mean profile','smoothed')
% 
%     figure()
%     new_os=mean_ap_os+mean_avk_os.*(mean_os-mean_ap_os);
%     plot(mean(lut_prof_os(:,3:end)),alt,'b-'), hold on
%     plot(mean(osiris_no2.num_dens),12.5:31.5,'bx-')
%     plot(new_os,12.5:39.5,'kx-')
%     title('OS mean smoothed')
%     legend('mean a priori','mean profile','smoothed')

end

