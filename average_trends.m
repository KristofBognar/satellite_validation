function average_trends()
% average trends for satellite instruments (after Hubert et al., 2016)

data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
load(data_file)
bruker_no2.tot_col=bruker_no2.part_col;

dt=12;

version=1;

if version==1 % use variance weigthed mean and corresponding uncertainty (Hubert)
    ncols=1;
elseif version==2 % use arithmetic mean, and uncertainty from mean bootstrap distribution (Gardiner)
    ncols=1000;
end
if 0
    
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


% brewer_o3_ds(brewer_o3_ds.sza>76.3,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OS_o3_trend=NaN(1,5);
OS_o3_trend_sig=NaN(5,ncols);

AF_o3_trend=NaN(1,4);
AF_o3_trend_sig=NaN(4,ncols);

AM_o3_trend=NaN(1,4);
AM_o3_trend_sig=NaN(4,ncols);

OS_no2_trend=NaN(1,4);
OS_no2_trend_sig=NaN(4,ncols);

AF_no2_trend=NaN(1,4);
AF_no2_trend_sig=NaN(4,ncols);

%%%

[OS_o3_trend(1),OS_o3_trend_sig(1,:)]=trend_time(osiris_o3_bw, brewer_o3_ds, version);
[OS_o3_trend(2),OS_o3_trend_sig(2,:)]=trend_time(osiris_o3_bk, bruker_o3, version);
[OS_o3_trend(3),OS_o3_trend_sig(3,:)]=trend_time(osiris_o3_bk, paris_o3, version);
[OS_o3_trend(4),OS_o3_trend_sig(4,:)]=trend_time(osiris_o3, gbs_o3, version);
[OS_o3_trend(5),OS_o3_trend_sig(5,:)]=trend_time(osiris_o3, saoz_o3, version);

[AF_o3_trend(1),AF_o3_trend_sig(1,:)]=trend_time(ace_fts_o3_bk, bruker_o3, version);
[AF_o3_trend(2),AF_o3_trend_sig(2,:)]=trend_time(ace_fts_o3_bk, paris_o3, version);
[AF_o3_trend(3),AF_o3_trend_sig(3,:)]=trend_twilight(ace_fts_o3, gbs_o3, version);
[AF_o3_trend(4),AF_o3_trend_sig(4,:)]=trend_twilight(ace_fts_o3, saoz_o3, version);

[AM_o3_trend(1),AM_o3_trend_sig(1,:)]=trend_time(ace_mae_o3_bk, bruker_o3, version);
[AM_o3_trend(2),AM_o3_trend_sig(2,:)]=trend_time(ace_mae_o3_bk, paris_o3, version);
[AM_o3_trend(3),AM_o3_trend_sig(3,:)]=trend_twilight(ace_mae_o3, gbs_o3, version);
[AM_o3_trend(4),AM_o3_trend_sig(4,:)]=trend_twilight(ace_mae_o3, saoz_o3, version);

%%%
[OS_no2_trend(1),OS_no2_trend_sig(1,:)]=trend_time(osiris_no2_bk, bruker_no2, version);
[OS_no2_trend(2),OS_no2_trend_sig(2,:)]=trend_time(osiris_no2, saoz_no2, version);
[OS_no2_trend(3),OS_no2_trend_sig(3,:)]=trend_time(osiris_no2, gbs_no2, version);
[OS_no2_trend(4),OS_no2_trend_sig(4,:)]=trend_time(osiris_no2uv, gbs_no2uv, version);

[AF_no2_trend(1),AF_no2_trend_sig(1,:)]=trend_time(ace_fts_no2_bk, bruker_no2, version);
[AF_no2_trend(2),AF_no2_trend_sig(2,:)]=trend_twilight(ace_fts_no2, saoz_no2, version);
[AF_no2_trend(3),AF_no2_trend_sig(3,:)]=trend_twilight(ace_fts_no2, gbs_no2, version);
[AF_no2_trend(4),AF_no2_trend_sig(4,:)]=trend_twilight(ace_fts_no2uv, gbs_no2uv, version);

if version==1
    
    OS_o3_trend_sig=OS_o3_trend_sig';
    AF_o3_trend_sig=AF_o3_trend_sig';
    AM_o3_trend_sig=AM_o3_trend_sig';

    OS_no2_trend_sig=OS_no2_trend_sig';
    AF_no2_trend_sig=AF_no2_trend_sig';
    
    
    %%% mean OSIRIS O3 trend

    n=length(OS_o3_trend);
    w=OS_o3_trend_sig.^(-2);

    OS_o3 = sum(OS_o3_trend .* w) / sum(w);
    OS_o3_sig = 1 / sqrt(sum(w));

    % OS_o3_sig = sqrt(  ( n/( (n-1)*sum(w)^2 ) ) * sum( w.^2 .* (OS_o3_trend-OS_o3).^2)  );
%     OS_o3_chi2=( sum( (OS_o3_trend-OS_o3).^2 ./ OS_o3_trend_sig.^2 ) ) / (n-1)

    %%% mean ACE-FTS O3 trend

    n=length(AF_o3_trend);
    w=AF_o3_trend_sig.^(-2);

    AF_o3 = sum(AF_o3_trend .* w) / sum(w);
    AF_o3_sig = 1 / sqrt(sum(w));

    % AF_o3_sig = sqrt(  ( n/( (n-1)*sum(w)^2 ) ) * sum( w.^2 .* (AF_o3_trend-AF_o3).^2)  );
%     AF_o3_chi2=( sum( (AF_o3_trend-AF_o3).^2 ./ AF_o3_trend_sig.^2 ) ) / (n-1)

    %%% mean ACE-MAESTRO O3 trend

    n=length(AM_o3_trend);
    w=AM_o3_trend_sig.^(-2);

    AM_o3 = sum(AM_o3_trend .* w) / sum(w);
    AM_o3_sig = 1 / sqrt(sum(w));

    % AM_o3_sig = sqrt(  ( n/( (n-1)*sum(w)^2 ) ) * sum( w.^2 .* (AM_o3_trend-AM_o3).^2)  );
%     AM_o3_chi2=( sum( (AM_o3_trend-AM_o3).^2 ./ AM_o3_trend_sig.^2 ) ) / (n-1)

    %%% mean OSIRIS NO2 trend

    n=length(OS_no2_trend);
    w=OS_no2_trend_sig.^(-2);

    OS_no2 = sum(OS_no2_trend .*w ) / sum(w);
    OS_no2_sig = 1 / sqrt(sum(w));

    % OS_no2_sig = sqrt(  ( n/( (n-1)*sum(w)^2 ) ) * sum( w.^2 .* (OS_no2_trend-OS_no2).^2)  );
%     OS_no2_chi2=( sum( (OS_no2_trend-OS_no2).^2 ./ OS_no2_trend_sig.^2 ) ) / (n-1)

    %%% mean ACE-FTS NO2 trend

    n=length(AF_no2_trend);
    w=AF_no2_trend_sig.^(-2);

    AF_no2 = sum(AF_no2_trend .* w) / sum(w);
    AF_no2_sig = 1 / sqrt(sum(w));

    % AF_no2_sig = sqrt(  ( n/( (n-1)*sum(w)^2 ) ) * sum( w.^2 .* (AF_no2_trend-AF_no2).^2)  );
%     AF_no2_chi2=( sum( (AF_no2_trend-AF_no2).^2 ./ AF_no2_trend_sig.^2 ) ) / (n-1)

elseif version==2
    
    OS_o3=mean(OS_o3_trend);
    OS_o3_sig=std(sum(OS_o3_trend_sig,1))*2;
    
    AF_o3=mean(AF_o3_trend);
    AF_o3_sig=std(sum(AF_o3_trend_sig,1))*2;
 
    AM_o3=mean(AM_o3_trend);
    AM_o3_sig=std(sum(AM_o3_trend_sig,1))*2;

    OS_no2=mean(OS_no2_trend);
    OS_no2_sig=std(sum(OS_no2_trend_sig,1))*2;
    
    AF_no2=mean(AF_no2_trend);
    AF_no2_sig=std(sum(AF_no2_trend_sig,1))*2;
    
end

disp(sprintf('OSIRIS O3 trend: %.1f +- %.1f %%/decade', [OS_o3, OS_o3_sig]));
disp(sprintf('ACE-FTS O3 trend: %.1f +- %.1f %%/decade', [AF_o3, AF_o3_sig]));
disp(sprintf('ACE-MAESTRO O3 trend: %.1f +- %.1f %%/decade', [AM_o3, AM_o3_sig]));

disp(sprintf('OSIRIS NO2 trend: %.1f +- %.1f %%/decade', [OS_no2, OS_no2_sig]));
disp(sprintf('ACE-FTS NO2 trend: %.1f +- %.1f %%/decade', [AF_no2, AF_no2_sig]));
    


end

function [trend, trend_sig] = trend_time(sat, gb, ver)

    dt=12;
    
    [ x1, x2, ~, ~, times] = find_coincidences_time( sat, gb, dt, false );
    
    diff=rel_abs_diff('rel',x1,x2);
    
    if max(diff)>500
        diff=diff./1e14;
    end

    if ver==1
        
        rfit=TrendAnalysis(times,diff);

        trend=rfit.trend*10;
        trend_sig=rfit.trend_sig*rfit.corr_factor*10;

    elseif ver==2
        
        rfit=TrendAnalysis(times,diff,true);

        trend=rfit.trend*10;
        trend_sig=rfit.trend_dist';

    end

end

function [trend, trend_sig] = trend_twilight(sat, gb, ver)

    [ x1, x2, ~, ~, times] = find_coincidences_twilight( sat, gb, false );
    
    diff=rel_abs_diff('rel',x1,x2);
    
    if max(diff)>500
        diff=diff./1e14;
    end

    if ver==1
        
        rfit=TrendAnalysis(times,diff);

        trend=rfit.trend*10;
        trend_sig=rfit.trend_sig*rfit.corr_factor*10;

    elseif ver==2
        
        rfit=TrendAnalysis(times,diff,true);

        trend=rfit.trend*10;
        trend_sig=rfit.trend_dist';
        
    end
    
end