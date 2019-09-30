function [] = corr_plot( x1e1, x2e2, instr1, instr2, loc_str, partcol,...
                         title_str,text_size,label_size )
%corr_plot(x1e1,x2e2,instr1, instr2, loc_str) correlation plot with labels and values
%
%

% figure()

if nargin==6
    title_str='';
    text_size=8; % for default figures: 8
    label_size=10; % matlab default is 10
end

% % % change frame identifiers to upper case
% % title_str=upper(title_str);

% for presentations
% text_size=12; 
% label_size=15; 

%% units/limits based on magnitude
if max(x1e1(:,1))<1000 % o3 in DU
    unitstr=' (DU)';

    if ~partcol
        limit=[200,650];
        xt=[200,400,600];
    else
        limit=[150,450];
        xt=[200,300,400];
    end
    
else % NO2 in 10^15 molec/cm2
    unitstr=' (x10^1^5 molec/cm^2)';
    
    limit=[0,6.5];
    xt=[0,2,4,6];
    
    % scale values
    mult=1e-15;
    x1e1=x1e1*mult;
    x2e2=x2e2*mult;
    
end

%% fit data

% Ls fit, no errors
[slope, y_int, R2, y_fit] = line_fit(x1e1(:,1), x2e2(:,1));
% % % [slope, y_int, R2, y_fit] = linear_regression(x1e1(:,1), x2e2(:,1), 0.31731, 0);

% fit with uniform weights, and with x,y errors
[y_int2(1), slope2(1), y_int2(2), slope2(2), R2_york]=...
    york_fit(x1e1(:,1)', x2e2(:,1)', 1, 1);

[y_int2e(1), slope2e(1), y_int2e(2), slope2e(2)]=...
    york_fit(x1e1(:,1)', x2e2(:,1)',x1e1(:,2)', x2e2(:,2)');

%% plot results

% 1 to 1 line
plot([limit],[limit],'k-','linewidth',1.4), hold on

% data points
% plot(x1e1(:,1),x2e2(:,1),'ro','markersize',5,'linewidth',1.1)
plot(x1e1(:,1),x2e2(:,1),'o','color',[0 0.6 0.6],'markersize',5,'linewidth',1.1)

% y= a*x + b, no errors
plot([limit], [slope(1)*limit(1)+y_int(1),slope(1)*limit(2)+y_int(1)],...
    'r--','linewidth',1.4)

% y= a*x + b, york fit with uniform errors
plot([limit], [slope2(1)*limit(1)+y_int2(1),slope2(1)*limit(2)+y_int2(1)],...
    'b--','linewidth',1.4)

% % % y= a*x + b, york fit with actual errors
% % plot([limit], [slope2e(1)*limit(1)+y_int2e(1),slope2e(1)*limit(2)+y_int2e(1)],...
% %     'g--','linewidth',1.4)

% force square plots
pbaspect([1 1 1])

% plot limits    
xlim([limit])
ylim([limit])

% uniform axes ticks
set(gca,'XTick',xt)
set(gca,'YTick',xt)

%% create labels

if nargin==2
    return
elseif nargin==4
    loc_str='bottomleft';
end

switch loc_str
    case 'left'
        ylabel([instr2 unitstr])
    case 'bottom'
        xlabel([instr1 unitstr])
    case 'bottomleft'
        ylabel([instr2 unitstr])
        xlabel([instr1 unitstr])
end
   
title(title_str,'FontWeight','Bold')

set(gca, 'FontSize', label_size)


% pring info on figure

% % text(0.5,0.1,['R^2=' num2str(round(R2,2))],'color','k','Units','normalized',...
% %     'fontsize',text_size)

text(0.05,0.93,sprintf('m=%.3f, b=%.2g', [slope(1),y_int(1)]),...
    'color','r','Units','normalized','fontsize',text_size)

text(0.05,0.83,sprintf('m=%.3f, b=%.2g', [slope2(1),y_int2(1)]),...
    'color','b','Units','normalized','fontsize',text_size)

% % text(0.05,0.85,['m=' num2str(round(slope2e(1),3)) ', b='...
% %     num2str(y_int2e(1)/mult,2)], 'color','g','Units','normalized',...
% %     'fontsize',8)

text(0.95,0.17,sprintf('N=%.0f', size(x1e1,1)), 'color','k','Units','normalized',...
    'fontsize',text_size,'HorizontalAlignment','right')

text(0.95,0.07,sprintf('{\\itR}=%.2f', sqrt(R2)),'color','k','Units','normalized',...
    'fontsize',text_size,'HorizontalAlignment','right')



%%% print error on fit parameters
% % text(0.05,0.9,{sprintf(['m=%.2f ' 177 ' %.2f'], [slope(1),slope(2)]),...
% %     sprintf(['b=%.2f ' 177 ' %.2f' ],[y_int(1),y_int(2)])},...
% %     'color','r','Units','normalized','fontsize',text_size)
% % 
% % text(0.05,0.72,{sprintf(['m=%.2f ' 177 ' %.2f'], [slope2(1),slope2(2)]),...
% %     sprintf(['b=%.2f ' 177 ' %.2f' ],[y_int2(1),y_int2(2)])},...
% %     'color','b','Units','normalized','fontsize',text_size)
% % 
% % text(0.95,0.17,['N=' num2str(size(x1e1,1))], 'color','k','Units','normalized',...
% %     'fontsize',text_size,'HorizontalAlignment','right')
% % 
% % text(0.95,0.07,sprintf('R=%.2f', sqrt(R2)),...
% %     'color','k','Units','normalized','fontsize',text_size,...
% %     'HorizontalAlignment','right')


end

