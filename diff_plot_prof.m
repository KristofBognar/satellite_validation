function diff_plot_prof( ax_in, prof1,prof2, instr1,instr2,...
                         loc_str,title_str,text_size,label_size,double_plot,alt)
%DIFF_PLOT_PROF plot absolute or relative differences between profiles prof1 and prof2

% altitude range 
if nargin<11, alt=14.5:51.5; end

% plot line width (dashed lines still 1)
lw=1.5;

prof1=prof1*1e-12;
prof2=prof2*1e-12;

% get differences
for i=1:length(alt)
    [abs_diff(i),abs_diff_sig(i)]=mean_diff('abs',prof1(:,i),prof2(:,i));
    [rel_diff(i),rel_diff_sig(i)]=mean_diff('rel',prof1(:,i),prof2(:,i));
end

disp(['Abs. diff.: ' num2str(mean(abs_diff)) 'e+12'])
disp(['Rel. diff.: ' num2str(mean(rel_diff)) ' %'])

if ~double_plot
    plot_color1=[0.6 0.6 0];
    plot_color2=[0.6 0 0.6];
    plot_color3=[0 0.6 0.6];
end


if ~double_plot

    axes(ax_in(1))
    
    
    plot(mean(prof1), alt, 'color', plot_color1,'linewidth',lw), hold on
    plot(mean(prof2), alt, 'color', plot_color2,'linewidth',lw)
    
    lgnd=legend(instr1,instr2);
    set(lgnd,'color','none');
    legend boxoff
    
%     text(0.87,0.7,sprintf('N=%.0f', size(prof1,1)),...
%          'color','k','Units','normalized','HorizontalAlignment','right')


    plot(mean(prof1)+std(prof1), alt, '--', 'color', plot_color1)
    plot(mean(prof1)-std(prof1), alt, '--', 'color', plot_color1)

    plot(mean(prof2)+std(prof2), alt, '--', 'color', plot_color2)
    plot(mean(prof2)-std(prof2), alt, '--', 'color', plot_color2)

    if strcmp(loc_str,'bottomleft')
        ylabel('Altitude (km)')
        xlabel({'Ozone profile','(x10^{12} molec/cm^3)'})
    elseif strcmp(loc_str,'left')
        ylabel('Altitude (km)')
    elseif strcmp(loc_str,'bottom')
        xlabel('Conc. (molec/cm^2)')
    end
    
    xlim([0,8])
    ylim([10,55])
    
    set(gca,'XTick',[0,2,4,6])
    set(gca,'TickLength',[0.02,0.025]);    
    
    set(gca,'XGrid','on')
    
    box off
    
    set(gca, 'FontSize', label_size)


    % abs diffs plot
    axes(ax_in(2))

    plot([0,0], [0,100], 'k-'), hold on

    plot(abs_diff,alt,'-','color',plot_color3,'linewidth',lw)
    plot(abs_diff+abs_diff_sig,alt,'--','color',plot_color3)
    plot(abs_diff-abs_diff_sig,alt,'--','color',plot_color3)

    if strcmp(loc_str,'bottomleft') || strcmp(loc_str,'bottom')
        xlabel({'\Delta_{abs}','(x10^{12} molec/cm^3)'})
    end

%     xlim([-1.4,1.4])
    xlim([-2.3,2.3])
    ylim([10,55])
    
    set(gca,'XGrid','on')
    set(gca,'YTickLabels',{})
    set(gca,'XTick',[-2,-1,0,1,2])
    set(gca,'TickLength',[0.02,0.025]);    
    
    box off
    
    text(0.6,0.88,[instr1 ' - ' instr2], 'color','k','Units','normalized',...
        'fontsize',text_size,'fontweight','normal')

    title(title_str,'FontWeight','Bold')

    set(gca, 'FontSize', label_size)
    
    
    % rel diffs plot
    axes(ax_in(3))

    plot([0,0], [0,100], 'k-'), hold on

    plot(rel_diff,alt,'-','color',plot_color3,'linewidth',lw)
    plot(rel_diff+rel_diff_sig,alt,'--','color',plot_color3)
    plot(rel_diff-rel_diff_sig,alt,'--','color',plot_color3)
    
    if strcmp(loc_str,'bottomleft') || strcmp(loc_str,'bottom')
        xlabel('\Delta_{rel} (%)')
    end
    
    xlim([-48,48])
    ylim([10,55])
    
    set(gca,'XGrid','on')
    set(gca,'YTickLabels',{})
    set(gca,'TickLength',[0.02,0.025]);
    
    box off
    set(gca, 'FontSize', label_size)    


% 
%     text(0.1,0.05,['\Delta_{abs}: '...
%         sprintf(['%3.3g ' 177 ' %3.3g\n'], [mean_abs_diff,mean_abs_diff_sig])],...
%         'color','k','Units','normalized','fontsize',text_size)
% 
%     text(0.9,0.05,['\Delta_{rel}: '...
%         sprintf(['%.1f ' 177 ' %.1f %%\n'], [mean_rel_diff,mean_rel_diff_sig])],...
%         'color','k','Units','normalized','fontsize',text_size,...
%         'HorizontalAlignment','right')
end    











end

