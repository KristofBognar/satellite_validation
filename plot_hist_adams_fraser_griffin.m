function plot_hist_adams_fraser_griffin(tg)

% ozone
% lists of differences (OS-GV, OS-SA, OS-BK, OS-PA, OS-BW,...
%                       AF-GV, AF-SA, AF-BK, AF-PA, AM-GV, AM-SA, AM-BK, AM-PA)
o3_current=[4.4,2.3,-2.1,-0.1,2.7,2.6,-0.5,-7.5,-4.3,-1.2,-4.4,-12,-8.8];
o3_current_err=[0.1,0.1,0.1,0.1,0.04,0.3,0.3,0.2,0.1,0.5,0.4,0.3,0.1];

o3_adams=[5.7,7.3,0.1,NaN,2.8,6.5,4.8,-4.7,NaN,5.0,1.6,-6.1,NaN];
o3_adams_err=[0.1,0.1,0.1,NaN,0.1,0.5,0.5,0.3,NaN,0.7,1,0.6,NaN];

% no2
% lists of differences (OS-GV, OS-GU, OS-SA, OS-BK, AF-GV, AF-GU, AF-SA, AF-BK)
no2_current=[-19.9,-8.1,-11.3,5.5,15.3,8.8,0.7,33.2];
no2_current_err=[0.3,0.4,0.6,0.3,0.9,1.1,1.1,0.6];
no2_adams=[-7.8,-3.3,10.2,12.2,15.2,13.6,12.7,NaN];
no2_adams_err=[0.4,0.3,0.9,0.2,1.3,2.2,2.1,NaN];

% sat-sat o3 (OS-AF, OS-AM, AF-AM)
sat_current=[1.2,6.7,5.9];
sat_current_err=[0.2,0.3,0.3];

sat_adams=[1.2,2.8,2.8];
sat_adams_err=[0.2,0.4,0.6];

if tg==99
    %% subplots

    figure()
    ax1=subplot(311);
    
    b1 = bar([o3_current', o3_adams'],0.7); hold on

%     b1(1).FaceColor = [0.6 0 0.6];
%     b1(2).FaceColor = [0.6 0.6 0];
    b1(1).FaceColor = [0.4 0.2 0.8];
    b1(2).FaceColor = [1 1 0];


    c = {'OS-GV','OS-SA','OS-BK','OS-PA','OS-BW','AF-GV','AF-SA','AF-BK','AF-PA',...
         'AM-GV','AM-SA','AM-BK','AM-PA'};

    set(gca,'xticklabel',c)

    % add other results

    % Fraser et al 2008
    plot([5.75,6.25],[3.2,3.2],'color',[0 .6 .6],'linewidth',2)
    plot([7.75,8.25],[-5.6,-5.6],'color',[.7 0 .7],'linewidth',2) % squeeze in Batchelor et al for labeling
    plot([8.75,9.25],[-5.2,-5.2],'color',[.25 0 .25],'linewidth',2) % squeeze in Fu et al AF-PA
    plot([7.75,8.25],[-3.6,-3.6],'color',[.6 .6 0],'linewidth',2) % squeeze in Griffin et al for labeling
    plot([8.75,9.25],[-3.5,-3.5],'color',[.6 .6 0],'linewidth',2) % squeeze in Griffin et al AF-PA

    % Fraser et al 2008, cont'd
    plot([5.75,6.25],[6.3,6.3],'color',[0 .6 .6],'linewidth',2)

    plot([6.75,7.25],[0.1,0.1],'color',[0 .6 .6],'linewidth',2)
    plot([6.75,7.25],[4.3,4.3],'color',[0 .6 .6],'linewidth',2)

    % plot([7.75,8.25],[-19.4,-19.4],'r-','linewidth',2)
    plot([9.75,10.25],[-1.2,-1.2],'color',[0 .6 .6],'linewidth',2)
    plot([9.75,10.25],[-2.6,-2.6],'color',[0 .6 .6],'linewidth',2)

    quiver(10,-11,0,-1,0,'linewidth',2,'maxheadsize',0.7,'color',[0 .6 .6])
    text(10,-10.5,'-19.4','horizontalalignment','center','fontsize',9)


    plot([10.75,11.25],[-12.9,-12.9],'color',[0 .6 .6],'linewidth',2)
    plot([10.75,11.25],[-1.9,-1.9],'color',[0 .6 .6],'linewidth',2)


    ll=legend('Current','Adams et al. (2012)', 'Fraser et al. (2008)',...
              'Batchelor et al. (2010)', 'Fu et al. (2011)','Griffin et al. (2017)',...
              'location','southwest');

    % add error bars (none for Fraser); need stupid capsize workaround...
    hErr=errorbar([1:13]-0.144,o3_current,o3_current_err, 'k.','marker','none');
    capsize(hErr,[1:13]-0.144,0.55);
    
%     hErr=errorbar([1:13]+0.144,o3_adams,o3_adams_err, 'k.','marker','none');
    ind=[1,2,3,5,6,7,8,10,11,12];
    hErr=errorbar(ind+0.144,o3_adams(ind),o3_adams_err(ind), 'k.','marker','none');
    capsize(hErr,ind+0.144,0.6);
    
    errorbar(8,-5.6,1.1,'color',[.6 0 .6],'marker','none'); % batchelor
    errorbar(9,-5.2,2.1,'color',[.3 0 .3],'marker','none'); % fu
    errorbar(8,-3.6,0.6,'color',[.6 .6 0],'marker','none'); % griffin
    errorbar(9,-3.5,0.6,'color',[.6 .6 0],'marker','none'); % griffin

    plot([5.5,5.5],[-15,10],'k-')
    plot([9.5,9.5],[-15,10],'k-')
    
    text(0.7,7,'a)','fontweight','bold');

%     xlabel('Instrument pairs')
    ylabel('Ozone \Delta_{rel} (%)')

    xlim([0.5,13.5])
    ylim([-13.1,8])

    set(gca,'XTickLabelRotation',45,'YGrid','on')
    set(gca, 'FontSize', 11)
    
    %%%%%%%%%%%%%%%%% NO2 %%%%%%%%%%%%%%%%%%%
    
    ax2=subplot(312);
     
    b2 = bar([no2_current', no2_adams'],0.7); hold on

    b2(1).FaceColor = [0.4 0.2 0.8];
    b2(2).FaceColor = [1 1 0];
    
    c = {'OS-GV','OS-GU','OS-SA','OS-BK','AF-GV','AF-GU','AF-SA','AF-BK'};

    set(gca,'xticklabel',c,'TickLength',[0.014,0.01])

    % Fraser et al 2008
    plot([4.75,5.25],[-10.7,-10.7],'color',[0 .6 .6],'linewidth',2)
    plot([4.75,5.25],[-13.9,-13.9],'color',[0 .6 .6],'linewidth',2)
    plot([4.75,5.25],[-19.7,-19.7],'color',[0 .6 .6],'linewidth',2)

    plot([6.75,7.25],[-11.9,-11.9],'color',[0 .6 .6],'linewidth',2)
    plot([6.75,7.25],[-13.6,-13.6],'color',[0 .6 .6],'linewidth',2)

%     legend('Current','Adams et al., 2012','Fraser et al., 2008','location','northwest')

    errorbar([1:8]-0.144,no2_current,no2_current_err, 'k.','marker','none');
    errorbar([1:8]+0.144,no2_adams,no2_adams_err, 'k.','marker','none');

    plot([4.5,4.5],[-30,40],'k-')
    
    text(0.7,32,'b)','fontweight','bold');
    
    xlabel('Instrument pairs')
    ylabel('NO_2 \Delta_{rel} (%)')

    xlim([0.5,8.5])
    ylim([-21,35])

    set(gca,'XTickLabelRotation',45,'YGrid','on')
    set(gca, 'FontSize', 11)

    %%%%%%%%%%%%%%%%% Sat O3 %%%%%%%%%%%%%%%%%%%
    
    ax3=subplot(313);
    
    b3 = bar([sat_current', sat_adams'],0.7); hold on

    b3(1).FaceColor = [0.4 0.2 0.8];
    b3(2).FaceColor = [1 1 0];
    
    c = {'OS-AF','OS-AM','AF-AM'};

    set(gca,'xticklabel',c,'TickLength',[0.026,0.01])
    
    hErr=errorbar([1:3]-0.144,sat_current,sat_current_err, 'k.','marker','none');
    capsize(hErr,[1:3]-0.144,3.2);

    hErr=errorbar([1:3]+0.144,sat_adams,sat_adams_err, 'k.','marker','none');
    capsize(hErr,[1:3]+0.144,3.2);
    
    plot([2.75,3.25],[5.5,5.5],'color',[0 .6 .6],'linewidth',2)
    plot([2.75,3.25],[7.8,7.8],'color',[0 .6 .6],'linewidth',2)

    quiver(3,8.4,0,0.6,0,'linewidth',1.8,'maxheadsize',1,'color',[0 .6 .6])
    text(2.74,8.7,'22.5','horizontalalignment','center','fontsize',9)
    
    text(0.7,8.93,'c)','fontweight','bold');
    
    xlabel('Instrument pairs')
    ylabel('Ozone \Delta_{rel} (%)')

    xlim([0.5,3.5])
    ylim([-1,9.5])

    set(gca,'XTickLabelRotation',45,'YGrid','on')
    set(gca, 'FontSize', 11)
    
    %%%%%%%%%%%%%% positions %%%%%%%%%%%%%%%%
    
    set(ax1,'Position',[0.09 0.58 0.86 0.38])
    set(ax2,'Position',[0.09 0.12 0.53 0.38])
    set(ax3,'Position',[0.751 0.12 0.199 0.38])
%     set(ll,'position',[0.7 0.22 0.22 0.18],'box','off')
    
    
elseif tg==1

    figure()
    hb = bar([o3_current', o3_adams'],0.7); hold on

    % hb(1).FaceColor = [0.6 0 0.6];
    % hb(2).FaceColor = [0.6 0.6 0];


    c = {'OS-GV','OS-SA','OS-BK','OS-PA','OS-BW','AF-GV','AF-SA','AF-BK','AF-PA',...
         'AM-GV','AM-SA','AM-BK','AM-PA'};

    set(gca,'xticklabel',c)

    % add other results

    % Fraser et al 2008
    plot([5.75,6.25],[3.2,3.2],'color',[0 .6 .6],'linewidth',2)
    plot([7.75,8.25],[-5.6,-5.6],'color',[.6 0 .6],'linewidth',2) % squeeze in Batchelor et al for labeling
    plot([8.75,9.25],[-5.2,-5.2],'color',[.3 0 .3],'linewidth',2) % squeeze in Fu et al AF-PA
    plot([7.75,8.25],[-3.6,-3.6],'color',[.6 .6 0],'linewidth',2) % squeeze in Griffin et al for labeling
    plot([8.75,9.25],[-3.5,-3.5],'color',[.6 .6 0],'linewidth',2) % squeeze in Griffin et al AF-PA

    % Fraser et al 2008, cont'd
    plot([5.75,6.25],[6.3,6.3],'color',[0 .6 .6],'linewidth',2)

    plot([6.75,7.25],[0.1,0.1],'color',[0 .6 .6],'linewidth',2)
    plot([6.75,7.25],[4.3,4.3],'color',[0 .6 .6],'linewidth',2)

    % plot([7.75,8.25],[-19.4,-19.4],'r-','linewidth',2)
    plot([9.75,10.25],[-1.2,-1.2],'color',[0 .6 .6],'linewidth',2)
    plot([9.75,10.25],[-2.6,-2.6],'color',[0 .6 .6],'linewidth',2)
    quiver(10,-11,0,-1,0,'linewidth',2,'maxheadsize',0.7,'color',[0 .6 .6])
    text(10,-10.5,'-19.4','horizontalalignment','center','fontsize',9)


    plot([10.75,11.25],[-12.9,-12.9],'color',[0 .6 .6],'linewidth',2)
    plot([10.75,11.25],[-1.9,-1.9],'color',[0 .6 .6],'linewidth',2)


    legend('Current','Adams et al., 2012', 'Fraser et al., 2008', 'Batchelor et al., 2010',...
           'Fu et al., 2011','Griffin et al., 2017','location','southwest')
    % legend('Current','Adams et al., 2012','location','southwest')


    % add error bars (none for Fraser)
    errorbar([1:13]-0.144,o3_current,o3_current_err, 'k.','marker','none');
    errorbar([1:13]+0.144,o3_adams,o3_adams_err, 'k.','marker','none');
    errorbar(8,-5.6,1.1,'color',[.6 0 .6],'marker','none') % batchelor
    errorbar(9,-5.2,2.1,'color',[.3 0 .3],'marker','none') % fu
    errorbar(8,-3.6,0.6,'color',[.6 .6 0],'marker','none') % griffin
    errorbar(9,-3.5,0.6,'color',[.6 .6 0],'marker','none') % griffin


    xlabel('Instrument pairs')
    ylabel('Ozone mean relative differences (%)')

    xlim([0.5,13.5])
    ylim([-13.1,8])

    set(gca,'XTickLabelRotation',45,'YGrid','on')


elseif tg==2
    figure()
    b = bar([no2_current', no2_adams'],0.7); hold on

    c = {'OS-GV','OS-GU','OS-SA','OS-BK','AF-GV','AF-GU','AF-SA','AF-BK'};

    set(gca,'xticklabel',c)

    % Fraser et al 2008
    plot([4.75,5.25],[-10.7,-10.7],'color',[0 .6 .6],'linewidth',2)
    plot([4.75,5.25],[-13.9,-13.9],'color',[0 .6 .6],'linewidth',2)
    plot([4.75,5.25],[-19.7,-19.7],'color',[0 .6 .6],'linewidth',2)

    plot([6.75,7.25],[-11.9,-11.9],'color',[0 .6 .6],'linewidth',2)
    plot([6.75,7.25],[-13.6,-13.6],'color',[0 .6 .6],'linewidth',2)

    legend('Current','Adams et al., 2012','Fraser et al., 2008','location','northwest')

    errorbar([1:8]-0.144,no2_current,no2_current_err, 'k.','marker','none');
    errorbar([1:8]+0.144,no2_adams,no2_adams_err, 'k.','marker','none');

    xlabel('Instrument pairs')
    ylabel('NO_2 mean relative differences (%)')

    xlim([0.5,8.5])
    ylim([-21,35])

    set(gca,'XTickLabelRotation',45,'YGrid','on')

end
end

function capsize(hErr,X,mult)
%     mult = 2;                               % twice as long
    b = hErr.Bar;                           % hidden property/handle
    drawnow                                 % populate b's properties
    vd = b.VertexData;
%     X=1:13;
    N = numel(X);                           % number of error bars
    capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
    newLength = capLength * mult;
    leftInds = N*2+1:2:N*6;
    rightInds = N*2+2:2:N*6;
    vd(1,leftInds,1) = [X-newLength, X-newLength];
    vd(1,rightInds,1) = [X+newLength, X+newLength];
    b.VertexData = vd;
end
