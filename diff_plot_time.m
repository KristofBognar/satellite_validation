function diff_plot_time( diff_type, x1,x2, param_in, plot_type, instr1,instr2,...
                         loc_str,title_str,text_size,label_size,double_plot)
%DIFF_PLOT plot absolute or relative differences between datasets x1 and x2

if nargin==7
    loc_str='';
    title_str='';
    text_size=8; % for default figures: 8
    label_size=10; % matlab default is 10
end
if nargin~=12
    double_plot=0;
end

% % % change frame identifiers to upper case
% % title_str=upper(title_str);

% get differences
diff=rel_abs_diff(diff_type,x1,x2);

% no2 or ozone?
o3_data=true;
if max(x1)>1000
    o3_data=false; 
    if strcmp(diff_type,'abs') diff=diff./1e14; end
    x1=x1./1e14;
    x2=x2./1e14;
end

% get mean diffs
[mean_abs_diff,mean_abs_diff_sig,mean_abs_diff_err]=mean_diff('abs',x1,x2);
[mean_rel_diff,mean_rel_diff_sig,mean_rel_diff_err]=mean_diff('rel',x1,x2);
rmsd=mean_diff('abs',x1,x2,'rmsd');

% if input is time, calculate ft and fractional year
if strcmp('timeseries',plot_type) || strcmp('seasonal',plot_type) ||...
   strcmp('ymean',plot_type) || strcmp('trends',plot_type)

    % get fractional time
    [ft,year]=mjd2k_to_ft(param_in);
    % get fractional year
    fracyear=year + ft./daysinyear(year);

end

% stuff for plotting multiple datasets
if o3_data
    title_x=0;
    title_y=[1.15,1.15];
    txt_offset=0.07;
    txt_x=[0.07,0.23,0.53,0.83];
else
    title_x=0;
    title_y=[1.15,1.15];
    txt_offset=0.07;
    txt_x=[0.07,0.23,0.53,0.83];
end

n_txt=0.26;

if ~double_plot
    plot_color=[0 0.6 0.6];

elseif double_plot==1 % first of two datasets
    plot_color=[0.6 0 0.6];
    
elseif double_plot==2 % second of two datasets
    plot_color=[0.6 0.6 0];
    txt_offset=-txt_offset;
    n_txt=1-n_txt;
    
elseif double_plot==11 % only one dataset
    plot_color=[0.6 0 0.6];
    txt_offset=0;

elseif double_plot==9 % only one dataset, o3 sat-sat plots
    plot_color=[0 0.6 0.6];
    title_y=[1.24,1.1];
    txt_offset=0;
    txt_x=[0.07,0,0.4,0.8];
    
elseif double_plot==12 % only one dataset, OS-BW plot
    plot_color=[0 0.6 0.6];
    txt_offset=0;
    txt_x(2)=0.25;
    
end

if any(double_plot==[9,12])
    text_color='k';
    text_weight='normal';
else
    text_color=plot_color;
    text_weight='bold';
end

switch plot_type
    case 'timeseries'

        plot([2003, 2018],[0,0],'k-'), hold on
        plot([2003, 2018],[mean(diff),mean(diff)],'--','color',plot_color), hold on
%         plot(fracyear, diff, 'r.')
        plot(fracyear, diff, '.','color',plot_color,'markersize',9)

%         [slope, y_int, R2] = line_fit(fracyear', diff);
%         fprintf('Slope: %3.3g +- %3.3g\n', slope(1),slope(2))

        xlim([2003,2018])
        
        if strcmp(loc_str,'bottom') || strcmp(loc_str,'bottomleft')
            xlabel('Year')
        end
        
        if o3_data
            ylim([-150,150])
        else
            ylim([-30,20])
        end
        
    case 'seasonal'

        plot([50, 300],[0,0],'k-'), hold on
        plot([50, 300],[mean(diff),mean(diff)],'--','color',plot_color), hold on
%         plot(fracyear, diff, 'r.')
        plot(ft+1, diff, '.','color',plot_color,'markersize',9)

%         plot(ft+1, diff, 'o','markerfacecolor',plot_color,'markersize',4,...
%             'markeredgecolor','w')

        xlim([50,300])
        
        if strcmp(loc_str,'bottom') || strcmp(loc_str,'bottomleft')
            xlabel('Day of the year')
        end
        
        if o3_data
            if strcmp(diff_type,'abs')
                ylim([-150,150])
            elseif strcmp(diff_type,'rel')
                ylim([-50,50])
            end
        else
            if strcmp(diff_type,'abs')
                ylim([-30,20])
            elseif strcmp(diff_type,'rel')
                ylim([-80,100])
            end
        end
        
    case 'ymean'
        % fit a least squares line to the yearly mean differences
        
        bins=unique(year); % get years with data
        yearplus=0.5; % for plotting mean values
        
        % set up arrays
        data=NaN(size(bins));
        err=NaN(size(bins));
        
        % loop over each year
        for i=1:length(bins)
        
            % find data for given year
            ind_year=find(year==bins(i));
%             ind_year=find(year==bins(i) & ft<90);
            
            % calculate mean and standard error for years with more than x
            % data points
            if length(ind_year)>10
                data(i)=nanmean(diff(ind_year));
%                 err(i)=nanstd(diff(ind_year)); % use std
                err(i)=nanstd(diff(ind_year))/sqrt(length(ind_year)); % use stdandard error
            end
        end
        
        % do least squares fit
        [slope, y_int] = line_fit(bins, data, err);
        
        % plot data
        plot(fracyear,diff,'.','color',[0,0,0]+0.75);  hold on
        
        plot(bins+yearplus, data, 's','color',plot_color,'linewidth',1.5),
        errorbar(bins+yearplus, data, err,'color',plot_color,'LineStyle','none');
        plot(bins+yearplus,slope(1)*(bins+yearplus)+y_int(1),'--','color',plot_color)
        
        trend=slope(1);
        trend_sig=slope(2);
        
        
% %         % do least squares fit to all data
% %         [slope, y_int] = line_fit(fracyear, diff);
% %         plot([2003, 2018],[0,0],'k-'), hold on
% %         plot(fracyear, diff, '.','color',plot_color,'markersize',9)
% %         plot(fracyear,slope(1)*fracyear+y_int(1),'--','color',plot_color)

    case 'trend'
        
        rfit=TrendAnalysis(param_in,diff);
        
        % plot data
        plot(rfit.x,rfit.y,'.','color',plot_color);  hold on
    
        % plot fit
        plot(rfit.x,rfit.trend*rfit.x+rfit.offset,'--','color',plot_color)
        
        trend=rfit.trend;
        trend_sig=rfit.trend_sig*rfit.corr_factor;

        trend_sig2=rfit.trend_sig;
%         trend_sig22=rfit.trend_sig_bootstrap;
        trend_sig3=rfit.trend_sig_weatherhead;
        
        
    case 'SZA'
        
        tmp_lim=min(param_in)-mod(min(param_in),5);

        plot([tmp_lim, 92],[0,0],'k-'), hold on
        plot([tmp_lim, 92],[mean(diff),mean(diff)],'--','color',plot_color), hold on

        plot(param_in, diff, '.','color',plot_color,'markersize',9)
        

        xlim([tmp_lim, 92])
        
        if strcmp(loc_str,'bottom') || strcmp(loc_str,'bottomleft')
            xlabel('SZA (degrees)')
        end
        
        if o3_data
            ylim([-150,150])
        else
            ylim([-30,20])
        end
        
        % adjust text so it's not covered by data points
        n_txt=0.07; % original: 0.26
        if double_plot==2, n_txt=1-n_txt; end
        
        % adjust more for OS-BW: single plot
        if double_plot==12
            title_x=-10; 
            n_txt=0.09;
            txt_x=txt_x-0.04;
        end
        
        

end
       
if strcmp(loc_str,'left') || strcmp(loc_str,'bottomleft')
    switch diff_type
        case 'abs'
            if o3_data
                ylabel('\Delta _{abs} (DU)')
            else
                ylabel('\Delta _{abs} (x10^{14} mol/cm^2)')
            end
        case 'rel'
            ylabel('\Delta_{rel} (%)')
    end
end



if ~double_plot
    title(title_str,'FontWeight','Bold')
    
    if ~o3_data, ylim([-30,30]); end
    
    set(gca, 'FontSize', label_size)

    text(0.1,0.9,[instr1 ' - ' instr2], 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','bold')

    text(0.1,0.05,['\Delta_{abs}: '...
        sprintf(['%.1f ' 177 ' %.1f\n'], [mean_abs_diff,mean_abs_diff_err])],...
        'color','k','Units','normalized','fontsize',text_size)

    text(0.9,0.05,['\Delta_{rel}: '...
        sprintf(['%.1f ' 177 ' %.1f %%\n'], [mean_rel_diff,mean_rel_diff_err])],...
        'color','k','Units','normalized','fontsize',text_size,...
        'HorizontalAlignment','right')
% %     text(0.9,0.05,['\Delta_{rel}: '...
% %         sprintf(['%.2f ' 177 ' %.2f %%\n'], [mean_rel_diff,mean_rel_diff_err])],...
% %         'color','k','Units','normalized','fontsize',text_size,...
% %         'HorizontalAlignment','right')
    
else
    
    box off
    
    text(title_x,title_y(1),title_str,'Units','normalized','FontWeight','Bold')
    
    set(gca, 'FontSize', label_size)
    set(gca, 'YGrid', 'on')
    
    text(txt_x(1),title_y(1)+txt_offset,[instr1 ' - ' instr2 ':'], 'color',text_color,...
        'Units','normalized','fontsize',text_size,'fontweight',text_weight)

    if strcmp('timeseries',plot_type) || strcmp('seasonal',plot_type)...
       || strcmp('SZA',plot_type)
        
    %         sprintf(['%3.3g ' 177 ' %3.3g\n'], [mean_abs_diff,mean_abs_diff_sig])],...
        text(txt_x(2),title_y(2)+txt_offset-0.06,['\Delta_{abs}: '...
            sprintf(['%.1f ' 177 ' %.1f\n'], [mean_abs_diff,mean_abs_diff_err])],...
            'color','k','Units','normalized','fontsize',text_size)

        text(txt_x(3),title_y(2)+txt_offset-0.06,['\Delta_{rel}: '...
            sprintf(['%.1f ' 177 ' %.1f %%\n'], [mean_rel_diff,mean_rel_diff_err])],...
            'color','k','Units','normalized','fontsize',text_size)
    % %     text(txt_x(3),title_y(2)+txt_offset-0.06,['\Delta_{rel}: '...
    % %         sprintf(['%.2f ' 177 ' %.2f %%\n'], [mean_rel_diff,mean_rel_diff_err])],...
    % %         'color','k','Units','normalized','fontsize',text_size)

        text(txt_x(4),title_y(2)+txt_offset,['{\itRMSD}: ' sprintf('%.1f', rmsd)],...
            'color','k','Units','normalized','fontsize',text_size)

        text(n_txt,0.07,sprintf('N=%.0f', length(diff)),'color',text_color,...
            'Units','normalized','fontsize',text_size,'fontweight',text_weight,...
            'HorizontalAlignment','center')
        
    else
        
        text(txt_x(2),title_y(2)+txt_offset-0.05,['m: '...
        sprintf('%.1f +- %.1f %%/decade\n',[trend*10,trend_sig*10])],...
        'color','k','Units','normalized','fontsize',text_size)

%         text(txt_x(2),title_y(2)+txt_offset-0.05,['m: '...
%         sprintf('%.1f +- %.1f, %.1f, %.1f, %.1f %%/decade\n',...
%         [trend,trend_sig,trend_sig2,trend_sig3]*10) ],...
%         'color','k','Units','normalized','fontsize',text_size)
    
        text(n_txt,0.07,sprintf('Need=%.0f yrs', ceil(rfit.years_needed)),...
            'color',text_color,...
            'Units','normalized','fontsize',text_size,'fontweight',text_weight,...
            'HorizontalAlignment','center')
    
    
    end

end    


end

