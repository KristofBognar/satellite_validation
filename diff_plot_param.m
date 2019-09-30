function diff_plot_param( diff_type, x1,x2, param, param_type, instr1, instr2,...
                          title_str,text_size,label_size )
%DIFF_PLOT plot absolute or relative differences between datasets x1 and x2
%as a function of solar zenith angle

if nargin==7
    title_str='';
    text_size=8; % for default figures: 8
    label_size=10; % matlab default is 10
end

% set plotting limits
plotlim=[floor(min(param)),ceil(max(param))];

% select widths for binning data
do_bin=true;
switch param_type
    case 'SZA'
        bin_width=1;
    case 'T'
        bin_width=5;
    case 'sPV'
        bin_width=0.2e-4;
        plotlim=[0,3e-4];
    case 'mean'
        do_bin=false;
        if max(param)<1000
            bin_width=20;
        else
            bin_width=0.2e15;
        end
    case 'day'
        bin_width=10;
        plotlim=[50,300];
    otherwise
        error('Parameter input not specified')
end

if isnan(plotlim(1)), return; end

% min number of counts required for a bing to be valid
% min_count=20;
min_count=5;


% get rel/abs difference of each pair
diff=rel_abs_diff(diff_type,x1,x2);

% get mean difference (same as mean(diff))
[mean_abs_diff,mean_abs_diff_sig]=mean_diff('abs',x1,x2);
[mean_rel_diff,mean_rel_diff_sig]=mean_diff('rel',x1,x2);


% bin data
if do_bin
    
    % lower limit: lowest multiple of bin_width
    llim=round( min(param) / bin_width ) * bin_width;
    % higher limit: highest multiple of bin_width
    hlim=round( max(param) / bin_width ) * bin_width;
    
    % list of bin centres
    bins=llim:bin_width:hlim; 
    
    % initialize binned data variables
    data=NaN(size(bins));
    err=NaN(size(bins));
    
    % loop over each bin
    for i=1:length(bins)
            
        % find coincidences that fall in bin
        ind_tmp=find(param>=bins(i)-bin_width/2 & param<(bins(i)+bin_width/2));
        
        % only include bin if there's sufficiant n.o. coincidences in it
        if length(ind_tmp)>min_count
            data(i)=nanmean(diff(ind_tmp));
            err(i)=nanstd(diff(ind_tmp));
        end
    end
    
    % plot unbinned data
    plot(param, diff,'.','color',[0,0.6,0.6]), hold on
    
    % plot mean and binned diff (incl. +-1 std error bars)
    plot(plotlim,[0,0],'k-'), hold on
    plot(plotlim,[mean(diff),mean(diff)],'k--'), hold on

    plot(bins, data, 's','color',[0 0.25 0.25])
    errorbar(bins, data, err,'color',[0 0.25 0.25], 'LineStyle','none');
    
end

xlabel(param_type)
if max(diff)<500
    ylabel(['\Delta_{' diff_type '} (DU)'])
else
    ylabel(['\Delta_{' diff_type '} (mol/cm^2)'])
end

title(title_str,'FontWeight','Bold')

set(gca, 'FontSize', label_size)


xlim(plotlim)

text(0.1,0.9,[instr1 ' - ' instr2], 'color','k','Units','normalized',...
    'fontsize',text_size,'fontweight','bold')

text(0.1,0.05,['\Delta_{abs}: ' sprintf('%3.3g\n', mean_abs_diff)],...
    'color','k','Units','normalized','fontsize',text_size)

text(0.9,0.05,['\Delta_{rel}: ' sprintf('%.1f %%\n', mean_rel_diff)],...
    'color','k','Units','normalized','fontsize',text_size,'HorizontalAlignment','right')


end

