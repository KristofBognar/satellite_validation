function [rfit_out, lfit_out, xfit_out] = TrendAnalysis(times, diff, do_bootstrap)
% a function to use both Weatherhead and Gardiner's emthod of anaylsis for
% trends in a dataset. The matrix inversion is a stronger fitting method
% under ideal cases
%
% Written by Paul Jeffery, modified by Kristof Bognar

%% Setup

if nargin==2, do_bootstrap=false; end

% unique days and indices of first occurrence of those days
[days,ind_days]=unique(floor(times));

% add last value to index array
ind_days=[ind_days; length(times)];

daily_diff=NaN(size(days));

% average daily values
for i=1:length(days), daily_diff(i)=mean(diff(ind_days(i):ind_days(i+1))); end

% Xdata=times;

% get fractional time
% [ft1,year1]=mjd2k_to_ft(times);
[ft2,year2]=mjd2k_to_ft(days);
% get fractional year
% fracyear_all=year1 + ft1./daysinyear(year1);
fracyear=year2 + ft2./daysinyear(year2);

%
Xdata=fracyear;

Ydata=daily_diff;

Seasonal='no'; % 'yes' for Gardiner et al. 2008, 'no' for Weatherhead et al. 1998

n_resamples=1000;

%% setup the model
if strcmp(Seasonal, 'yes') %as in gardiner
%     modelfit='u+w*x+b1*cos(2*pi*x)+b2*sin(2*pi*x)+b3*cos(2*pi*x)+b4*sin(2*pi*x)+b5*cos(2*pi*x)+b6*sin(2*pi*x)+b7*cos(2*pi*x)+b8*sin(2*pi*x)';
    modelfit='u+w*x+b1*cos(2*pi*x)+b2*sin(2*pi*x)+b3*cos(2*pi*x)+b4*sin(2*pi*x)+b5*cos(2*pi*x)+b6*sin(2*pi*x)';
elseif strcmp(Seasonal, 'no') %as in weatherhead
    modelfit='u+w*x+b1';
end
% no model for robust fit, 2-param line fit by default

%% fit the data

%%% using built-in fit (ordinary least squares)
warning off
lfit=fit(Xdata, Ydata, modelfit);

%%% using matrix inversion
%y in y=AX+e, where e is the model error (that we ignore), y the obs, A  the
%model and X the coeff (at+b0+b1cos2pit+b2sin2pit+etc

%generate the model matrix:
if strcmpi(Seasonal, 'yes')==1 %as in gardiner
    model=[Xdata ones(length(Xdata),1) cos(2*pi*Xdata) sin(2*pi*Xdata) ...
           cos(4*pi*Xdata) sin(4*pi*Xdata) cos(6*pi*Xdata) sin(6*pi*Xdata) ...
           cos(8*pi*Xdata) sin(8*pi*Xdata)];
else %as in weatherhead
    model=[Xdata ones(length(Xdata),1)];
end

model2=[Xdata ones(length(Xdata),1)];

%fit with mldivide syntax
xfit=model\Ydata;

%%% using built-in robust fit
% b + w*x, based on Holland and Welsch, 1977 (sigma estimate from Street, 1988)
[rfit,rparam]=robustfit(Xdata,Ydata);
rfit=flipud(rfit);


%% Bootstrap resampling
    
%calculate residuals R0/N
resid_lfit=Ydata-feval(lfit, Xdata);
resid_xfit=Ydata-model*xfit;
resid_rfit=Ydata-model2*rfit;

%generate storage arrays
if strcmpi(Seasonal, 'yes')==1 
    fitparam_lfit=nan(n_resamples, 10);
    fitparam_xfit=nan(n_resamples,10);
else
    fitparam_lfit=nan(n_resamples, 2);
    fitparam_xfit=nan(n_resamples,2);
end

fitparam_rfit=nan(n_resamples,2);

if do_bootstrap
    % loop over number of required resample fits
    n=0;
    for i=1:n_resamples 

        % display progress info
        disp_str=['Resample fit ', num2str(i), '/', num2str(n_resamples)];
        % stuff to delete last line and reprint updated message
        fprintf(repmat('\b',1,n));
        fprintf(disp_str);
        n=numel(disp_str);    

        % get random sample from residuals
        sample_ind=round(rand(length(Xdata),1)*(length(Xdata)-1))+1; 

        % calculate resampled datasets
        y_resample_lfit=feval(lfit, Xdata)+resid_lfit(sample_ind);
        y_resample_xfit=model*xfit+resid_xfit(sample_ind);
        y_resample_rfit=model2*rfit+resid_rfit(sample_ind);

        % do new fits
        lfit_tmp=fit(Xdata, y_resample_lfit, modelfit);
        xfit_tmp=model\y_resample_xfit;
        rfit_tmp=robustfit(Xdata,y_resample_rfit);

        % store least squares fit data
        if strcmpi(Seasonal, 'yes')==1
            fitparam_lfit(i,1)=lfit_tmp.w;
            fitparam_lfit(i,2)=lfit_tmp.u;
            fitparam_lfit(i,3)=lfit_tmp.b1;
            fitparam_lfit(i,4)=lfit_tmp.b2;
            fitparam_lfit(i,5)=lfit_tmp.b3;
            fitparam_lfit(i,6)=lfit_tmp.b4;
            fitparam_lfit(i,7)=lfit_tmp.b5;
            fitparam_lfit(i,8)=lfit_tmp.b6;
%             fitparam_lfit(i,9)=lfit_tmp.b7;
%             fitparam_lfit(i,10)=lfit_tmp.b8;
        else
            fitparam_lfit(i,1)=lfit_tmp.w;
            fitparam_lfit(i,2)=lfit_tmp.u+lfit_tmp.b1;
        end

        % store martix fit data
        fitparam_xfit(i,:)=xfit_tmp';

        % store robust fit data
        fitparam_rfit(i,:)=fliplr(rfit_tmp');

        clear sample_ind y_resample_lfit y_resample_xfit y_resample_rfit lfit_tmp xfit_tmp rfit_tmp

    end
    clear i
    
    % trend stats

    % calculate error on trend component using Bootstrap
    sw_lfitB=2*std(fitparam_lfit(:,1));
    sw_xfitB=2*std(fitparam_xfit(:,1));
    sw_rfitB=2*std(fitparam_rfit(:,1));

else

    sw_lfitB=NaN;
    sw_xfitB=NaN;
    sw_rfitB=NaN;
    
end


%% trend stats

% determine autocorrelation of the residuals (residuals assumed to represent the noise)
phi_lfit=corr(resid_lfit(2:end),resid_lfit(1:end-1));
phi_xfit=corr(resid_xfit(2:end),resid_xfit(1:end-1));
phi_rfit=corr(resid_rfit(2:end),resid_rfit(1:end-1));
   
% determine error on noise sn using method from weatherhead
sn_lfitW=std(resid_lfit);
sn_xfitW=std(resid_xfit);
sn_rfitW=std(resid_rfit);

% number of years to detect a trend using Weatherhead eqn 3
years_lfit=((3.3*sn_lfitW/(abs(lfit.w)))*sqrt((1+phi_lfit)/(1-phi_lfit)))^(2/3);
years_xfit=((3.3*sn_xfitW/abs(xfit(1)))*sqrt((1+phi_xfit)/(1-phi_xfit)))^(2/3);
years_rfit=((3.3*sn_rfitW/abs(rfit(1)))*sqrt((1+phi_rfit)/(1-phi_rfit)))^(2/3);

% calculate error on trend component using Weatherhead eqn 2
sw_lfitW=(sn_lfitW/(Xdata(end)-Xdata(1))^(3/2))*sqrt((1+phi_lfit)/(1-phi_lfit));
sw_xfitW=(sn_xfitW/(Xdata(end)-Xdata(1))^(3/2))*sqrt((1+phi_xfit)/(1-phi_xfit));
sw_rfitW=(sn_rfitW/(Xdata(end)-Xdata(1))^(3/2))*sqrt((1+phi_rfit)/(1-phi_rfit));

%% output 

% ols fit
lfit_out=struct();

lfit_out.trend=lfit(1);
lfit_out.trend_sig_bootstrap=sw_lfitB;
lfit_out.trend_sig_weatherhead=sw_lfitW;

lfit_out.corr_factor=sqrt((1+phi_lfit)/(1-phi_lfit));
lfit_out.trend_dist=fitparam_lfit(:,1);

lfit_out.years_needed=years_lfit;

lfit_out.offset=lfit(2);

% matrix fit
xfit_out=struct();

xfit_out.trend=xfit(1);
xfit_out.trend_sig_bootstrap=sw_xfitB;
xfit_out.trend_sig_weatherhead=sw_xfitW;

xfit_out.corr_factor=sqrt((1+phi_xfit)/(1-phi_xfit));
xfit_out.trend_dist=fitparam_xfit(:,1);

xfit_out.years_needed=years_xfit;

xfit_out.offset=xfit(2);

%
rfit_out=struct();

rfit_out.x=fracyear;
rfit_out.y=daily_diff;

rfit_out.trend=rfit(1);
rfit_out.trend_sig=rparam.se(2)*2; % take 2 sigma
rfit_out.trend_sig_bootstrap=sw_rfitB;
rfit_out.trend_sig_weatherhead=sw_rfitW;

rfit_out.corr_factor=sqrt((1+phi_rfit)/(1-phi_rfit));

rfit_out.trend_dist=fitparam_rfit(:,1);
rfit_out.robust_param=rparam;

rfit_out.years_needed=years_rfit;

rfit_out.offset=rfit(2);
rfit_out.offset_sig=rparam.se(1)*2;

if do_bootstrap, fprintf('\n'); end

end
