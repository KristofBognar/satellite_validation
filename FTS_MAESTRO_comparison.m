function FTS_MAESTRO_comparison()
%FTS_MAESTRO_COMPARISON compare FTS and MAESTRO columns as a function of
%action table numbers /versions

%% setup
global alphabet text_size label_size ii

text_size=10; % matlab default is 8
label_size=10; % matlab default is 10
alphabet = 'abcdefghijklmnopqrstuvwxyz';

% load processed data
data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
load(data_file)

alt=14.5:51.5;

% AF--AM comparisons, no maestro error filer
% abs_all=15.3;
% abs_all_err=0.7;
% rel_all=5.9;
% rel_all_err=0.3;
[abs_all,abs_all_err,rel_all,rel_all_err,N_all]=...
    diff_twilight(ace_mae_o3,ace_fts_o3);

abs_all=round(abs_all,1);
abs_all_err=round(abs_all_err,1);
rel_all=round(rel_all,1);
rel_all_err=round(rel_all_err,1);

% find all maestro action tables
tables=unique(ace_mae_o3{:,19:20},'rows');

tables=sortrows(tables,[2,1]);

% set up differences table -- not used ATM
diffs_table=table;
diffs_table.action_table='All';
diffs_table.abs=abs_all;
diffs_table.abs_err=abs_all_err;
diffs_table.rel=rel_all;
diffs_table.rel_err=rel_all_err;
diffs_table.N=550;

figure()
nrow=4;
ncol=6;
fig_ax = tight_subplot(nrow,ncol,[0,0],[0.1,0.05],[0.1,0.05]);

%% loop through individual action tables
for i=1:size(tables,1)

    % find matching MAESTRO measurements
    ind=find(ace_mae_o3.action_table==tables(i,1) & ace_mae_o3.A0B1==tables(i,2));
    ace_mae_tmp=ace_mae_o3(ind,:);
    
    % make action table name
    AB='A';
    if tables(i,2)==1, AB='B'; end
    label=[AB  num2str(tables(i,1))];
    
    % find/save AF diffs
    [mean_abs_diff,mean_abs_diff_err,mean_rel_diff,mean_rel_diff_err,N,prof1,prof2]=...
        diff_twilight(ace_mae_tmp,ace_fts_o3);
    
    mean_abs_diff=round(mean_abs_diff, 1);
    mean_abs_diff_err=round(mean_abs_diff_err, 1);
    mean_rel_diff=round(mean_rel_diff, 1);
    mean_rel_diff_err=round(mean_rel_diff_err, 1);
    
    tmp_table=table;
    tmp_table.action_table=label;
    tmp_table.abs=mean_abs_diff;
    tmp_table.abs_err=mean_abs_diff_err;
    tmp_table.rel=mean_rel_diff;
    tmp_table.rel_err=mean_rel_diff_err;
    tmp_table.N=N;

    diffs_table=[diffs_table; tmp_table];
    
    % plot mean profiles
    axes(fig_ax(i));
    plot(mean(prof2), alt), hold on
    plot(mean(prof1), alt, '--'),
    
    xlim([0,7]*1e12)
    set(fig_ax(i),'XTick',[0 2 4 6]*1e12);
    if ~any(i==[(nrow-1)*ncol+1:nrow*ncol]), set(fig_ax(i),'XTickLabel',[]); end
    
    ylim([10,55])
    set(fig_ax(i),'YTick',[10 20 30 40 50]);    
    if ~any(i==[1:ncol:nrow*ncol]), set(fig_ax(i),'YTickLabel',[]); end
    set(fig_ax(i),'TickLength',[.03 .03]);
    
    box off
    grid on

    fill([7,7,1.2,1.2]*1e12,[35,55,55,35],'w','edgecolor','none')

    % highlight abs/rel diff if significantly different from mean diff
    % (within standard error)
    if abs(mean_abs_diff-abs_all)>(mean_abs_diff_err+abs_all_err)
        fill([7,7,1.2,1.2]*1e12,[41.5,48,48,41.5],'r','edgecolor','none','facealpha',0.2)
    end
    if abs(mean_rel_diff-rel_all)>(mean_rel_diff_err+rel_all_err)
        fill([7,7,1.2,1.2]*1e12,[35,41.5,41.5,35],'r','edgecolor','none','facealpha',0.2)
    end
    
    % add diffs and legend
%     text(0.25,0.92,['- - AF  - AM   ' label], 'color','k','Units','normalized',...
%         'fontsize',text_size)
    text(0.5,0.92,label, 'color','k','Units','normalized',...
        'fontsize',text_size)
    text(0.03,0.08,['N = ' num2str(N)], 'color','k','Units','normalized',...
        'fontsize',text_size)
    text(0.25,0.77,['\Delta_{abs}= ' num2str(mean_abs_diff) 177 num2str(mean_abs_diff_err) ' DU'],...
        'color','k','Units','normalized','fontsize',text_size)
    text(0.25,0.62,['\Delta_{rel} = ' num2str(mean_rel_diff) 177 num2str(mean_rel_diff_err) ' %'],...
        'color','k','Units','normalized','fontsize',text_size)
    
    
end

% add title and some axis labels
axes(fig_ax(3));
title(['All profiles (no error filter):  \Delta_{abs} = ' num2str(abs_all) ' ' 177 ...
       ' ' num2str(abs_all_err) ' DU;  \Delta_{rel} = ' num2str(rel_all) ' ' 177 ...
       ' ' num2str(rel_all_err) ' %;  N = ' num2str(N_all) ';  lines:  - - FTS  ' ...
       8211 ' MAESTRO'])
axes(fig_ax(1+ncol));
ylabel('Altitude (km)')
axes(fig_ax(nrow*ncol-floor(ncol/2)));
xlabel('molec/cm^3')

set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(gcf, 'Position', [100, 100, 185*ncol, 800]);


abs_changes=abs(diffs_table.abs-abs_all);
abs_singificant=zeros(size(diffs_table.abs));
abs_singificant(abs_changes>abs_all_err+diffs_table.abs_err)=1;
diffs_table.abs_significant=abs_singificant;

rel_changes=abs(diffs_table.rel-rel_all);
rel_singificant=zeros(size(diffs_table.rel));
rel_singificant(rel_changes>rel_all_err+diffs_table.rel_err)=1;
diffs_table.rel_significant=rel_singificant;

figure()
subplot(2,4,[1 2 3])
plot([0,26],[diffs_table.abs(1)+diffs_table.abs_err(1),...
             diffs_table.abs(1)+diffs_table.abs_err(1)],'k--'), hold on
plot([0,26],[diffs_table.abs(1)-diffs_table.abs_err(1),...
             diffs_table.abs(1)-diffs_table.abs_err(1)],'k--'), hold on

errorbar(diffs_table.abs,diffs_table.abs_err,'ks')

set(gca,'xtick',[1:length(diffs_table.abs)])
set(gca,'xticklabel',cellstr(diffs_table.action_table)','XTickLabelRotation',90)
xlim([0,26])
ylabel('\Delta_{abs} (DU)')
grid on

subplot(2,4,[5 6 7])
plot([0,26],[diffs_table.rel(1)+diffs_table.rel_err(1),...
             diffs_table.rel(1)+diffs_table.rel_err(1)],'k--'), hold on
plot([0,26],[diffs_table.rel(1)-diffs_table.rel_err(1),...
             diffs_table.rel(1)-diffs_table.rel_err(1)],'k--'), hold on
errorbar(diffs_table.rel,diffs_table.rel_err,'ks')
set(gca,'xtick',[1:length(diffs_table.abs)])
set(gca,'xticklabel',cellstr(diffs_table.action_table)','XTickLabelRotation',90)
xlim([0,26])
ylabel('\Delta_{rel} (%)')
grid on



subplot(244)
histogram(diffs_table.abs,15), hold on
plot([diffs_table.abs(1),diffs_table.abs(1)],[0,4],'k-','linewidth',1.3)
xlabel('\Delta_{abs} (DU)')
ylabel('#')

subplot(248)
histogram(diffs_table.rel,15), hold on
plot([diffs_table.rel(1),diffs_table.rel(1)],[0,6],'k-','linewidth',1.3)
xlabel('\Delta_{rel} (%)')
ylabel('#')


set(findall(gcf,'-property','FontName'),'FontName','Arial')
set(gcf, 'Position', [100, 100, 1110, 400]);

end

function [mean_abs_diff,mean_abs_diff_err,mean_rel_diff,mean_rel_diff_err,N,prof1,prof2]=...
    diff_twilight(table1, table2)

    global alphabet text_size label_size ii
    
    [ x1, x2, ~, ~, ~, ~, prof1, prof2] =...
        find_coincidences_twilight( table2, table1, true );
    
    [mean_abs_diff,mean_abs_diff_sig,mean_abs_diff_err]=mean_diff('abs',x1,x2);
    [mean_rel_diff,mean_rel_diff_sig,mean_rel_diff_err]=mean_diff('rel',x1,x2);
    
    N=length(x1);
    
end
