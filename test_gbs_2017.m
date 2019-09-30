% test 2017 gbs ozone data -- affected by amp signal

comps=1;
rms_and_err=0;

% tg='O3';
tg='NO2';

msize=10;
load('/home/kristof/work/satellite_validation/all_data_nosmooth.mat')


if comps==1

    [ x1, x2, e1, e2, times,~, DMP ] =...
            find_coincidences_time( brewer_o3_ds, gbs_o3, 12, false );
%     [ x1, x2, e1, e2, times, DMP ] =...
%             find_coincidences_twilight( saoz_o3, gbs_o3, false );

    [ft,year]=mjd2k_to_ft(times);

    ind=find(year==2017 & ft<119);

    % plot([200,550],[200,550],'b-'), hold on
    % 
    % plot(x1,x2,'k.'),
    % plot(x1(ind),x2(ind),'rx')

    figure()
    set(gcf, 'Position', [100, 100, 1000, 600]);
    msize=12;

    for i=2004:2017
        ind=find(year==i);

        plot([40,280],[0,0],'k-'), hold on
        plot([70.785,70.785],[-150,150],'k-')
        plot([119,119],[-150,150],'k--')

        plot(ft,x1-x2,'.','color',[1 1 1]*0.5,'markersize',msize)
        plot(ft(ind),x1(ind)-x2(ind),'k.','markersize',msize)

        title(num2str(i))
        xlim([50,280])

        xlabel('Fractional time (UTC)')
        ylabel('O_3 \Delta_{abs} (X-GBS)')

        if i<2017
            waitforbuttonpress;
            clf
        end

    end
end

if rms_and_err==1
    
    % read all RMS values
    ft_all=[];
    rms_all=[];
    yr_all=[];
    
    for i=2004:2017
        
        load(['/home/kristof/work/GBS/VCD_results/UT-GBS_' tg '_VCD_' num2str(i) '.mat']);
        ft_all=[ft_all; dscd_S.fd];
        rms_all=[rms_all; dscd_S.rms];
        yr_all=[yr_all; dscd_S.year];
        
    end

    if strcmp(tg,'O3')
        tg_table=gbs_o3;
        ytop=60;
    elseif strcmp(tg,'NO2')
        tg_table=gbs_no2;
        ytop=2e15;
    end
    
    error=sqrt(tg_table.sigma_mean_vcd.^2 + tg_table.std_vcd.^2);

    figure()
    set(gcf, 'Position', [100, 100, 1000, 750]);
    
    
    for i=2004:2017
        
        % plot QDOAS RMS
        subplot(211)
        hold on

        plot(0,0,'.','color',[1 1 1]*0.5,'markersize',msize)
        plot(0,0,'b.','markersize',msize)

        plot([70.785,70.785],[0,6e-3],'k-')
        plot([119,119],[0,6e-3],'k--')
        
        ind=find(yr_all==i);
        plot(ft_all-1,rms_all,'.','color',[1 1 1]*0.5,'markersize',msize)
        
        plot(ft_all(ind)-1,rms_all(ind),'k.','markersize',msize)

        legend('2004-2017','Title year','Shutter out, diffuser on','Lamp signal gone')
        ylabel([tg ' RMS'])
        xlim([40,300])
        ylim([0,4.5e-3])
        title(num2str(i))
        
        % plot total VCD error
        subplot(212)
        hold on

        plot([70.785,70.785],[0,ytop],'k-')
        plot([119,119],[0,ytop],'k--')
        
        ind=find(gbs_o3.year==i);
        plot(tg_table.fractional_time,error,'.','color',[1 1 1]*0.5,'markersize',msize)

        plot(tg_table.fractional_time(ind),error(ind),'k.','markersize',msize);
        
        xlabel('Fractional time (UTC)')
        ylabel([tg ' VCD error'])
        xlim([40,300])
    
        if i<2017
            waitforbuttonpress;
            clf
        end
        
        
    end
     
    
    
end
    
    
    

    
    
    