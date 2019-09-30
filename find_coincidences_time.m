function [ vcd1, vcd2, error1, error2, times, time_diff, DMP, prof1, prof2] = ...
         find_coincidences_time( table1, table2, dt, partcol )
%[ vcd1, vcd2, error1, error2 ] = FIND_COINCIDENCES_TIME(table1, table2, dt)
%   Find coincident measurements based on time constraint 
%
%   Input: 
%       table1, table2: tables of data, only required column is 'mjd2k'
%       (modified julian date with 2000 Jan. 1, 00:00 = 0)
%       dt: time window for coincidences (in hours), default is 12h
%       partcol (optional): if true, output is partial columns instead of
%       total columns (requires satellite files)
%
%   Output:
%       vcd1, vcd2, error1, error2: pairs of VCD and error values, assigned
%       based on expected input data (GBS, Brewer, satellites, etc)
%       time_diff: time difference between value pairs, in hours
%       sza1, sza2: pairs of solar zenith angle values corresponding to
%           vcd1 and vcd2
%       
%   Note: function eats memory and processing power if data dimensions are
%       in the 10 000s
%
% Kristof Bognar, January 2018

% if 0
% table1=saoz;
% table2=brewer_ds;
% else
% table1=brewer_ds;
% table2=saoz;
% end
% 
% dt=0.5;
% partcol=false;


% input arguments
if nargin<4
    partcol=false;
    if nargin<3
        dt=12;
    end
end

% convert time difference to days
dt=dt/24;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% find temporal coincidences

% find differences between each element
% bsxfun uses less memory (and requires fewer lines of code) compared to repmat
diff=abs(bsxfun(@minus,table2.mjd2k,table1.mjd2k'));
% each column represents the difference of one table1 time to all the table2 times
% each row represents the difference of one table2 time to all the table1 times

% indices of minima in table2, corresponding to each table1 element
% (closest times to each table1 time)
[min_val12,min_ind12]=min(diff,[],1);

% indices of minima in table1, corresponding to each table2 element
% (closest times to each table2 time)
[min_val21,min_ind21]=min(diff,[],2);

% save all coincidences: we want each value compared to the nearest meas. in
% the other dataset, and do this MUTUALLY (datasets are treated as equal)
% coincidences array has columns:
%   table1 indices
%   table2 indices
%   time difference (in mjd, so it's days)

% all table1 indices, with corresponding closest table2 values
coincidences=[[1:length(min_ind12)]', min_ind12', min_val12'];

% all table2 indices, with corresponding closest table1 values
coincidences=[coincidences;...
              [min_ind21, [1:length(min_ind21)]', min_val21]];

% discard coincidences outside time window
coincidences(coincidences(:,3)>dt,:)=[];

% remove double-counted entries
coincidences=unique(coincidences,'rows');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assign output values

% figure out date type for each table, based on unique columns
% 1: satelite, 2: bruker, 3: brewer 4: GBS

if any(strcmp('dist',table1.Properties.VariableNames)==1)
    type1=1;
elseif any(strcmp('tot_col_err_sys',table1.Properties.VariableNames)==1)
    type1=2;
elseif any(strcmp('tot_col_err',table1.Properties.VariableNames)==1)
    type1=3;
elseif any(strcmp('mean_vcd',table1.Properties.VariableNames)==1)
    type1=4;
end

if any(strcmp('dist',table2.Properties.VariableNames)==1)
    type2=1;
elseif any(strcmp('tot_col_err_sys',table2.Properties.VariableNames)==1)
    type2=2;
elseif any(strcmp('tot_col_err',table2.Properties.VariableNames)==1)
    type2=3;
elseif any(strcmp('mean_vcd',table2.Properties.VariableNames)==1)
    type2=4;
end

% assign VCD and error values
% table 1
if type1==1
    if ~partcol
        vcd1=table1.tot_col(coincidences(:,1));
        error1=table1.tot_col_err(coincidences(:,1));
    else
        vcd1=table1.part_col(coincidences(:,1));
        error1=table1.part_col_err(coincidences(:,1));
    end        
elseif type1==2
    if ~partcol
        vcd1=table1.tot_col(coincidences(:,1));
        error1=sqrt(table1.tot_col_err_sys(coincidences(:,1)).^2 +...
                    table1.tot_col_err_rand(coincidences(:,1)).^2 );
    else
        vcd1=table1.part_col(coincidences(:,1));
        error1=sqrt(table1.tot_col_err_sys(coincidences(:,1)).^2 +...
                    table1.tot_col_err_rand(coincidences(:,1)).^2 );
    end
elseif type1==3
    vcd1=table1.tot_col(coincidences(:,1));
    error1=table1.tot_col_err(coincidences(:,1));
elseif type1==4
    vcd1=table1.mean_vcd(coincidences(:,1));
    error1=sqrt(table1.sigma_mean_vcd(coincidences(:,1)).^2 +...
                table1.std_vcd(coincidences(:,1)).^2 );
end

% table 2
if type2==1
    if ~partcol
        vcd2=table2.tot_col(coincidences(:,2));
        error2=table2.tot_col_err(coincidences(:,2));
    else
        vcd2=table2.part_col(coincidences(:,2));
        error2=table2.part_col_err(coincidences(:,2));
    end        
elseif type2==2
    if ~partcol
        vcd2=table2.tot_col(coincidences(:,2));
        error2=sqrt(table2.tot_col_err_sys(coincidences(:,2)).^2 +...
                    table2.tot_col_err_rand(coincidences(:,2)).^2 );
    else
        vcd2=table2.part_col(coincidences(:,2));
        error2=sqrt(table2.tot_col_err_sys(coincidences(:,2)).^2 +...
                    table2.tot_col_err_rand(coincidences(:,2)).^2 );

    end
elseif type2==3
    vcd2=table2.tot_col(coincidences(:,2));
    error2=table2.tot_col_err(coincidences(:,2));
elseif type2==4
    vcd2=table2.mean_vcd(coincidences(:,2));
    error2=sqrt(table2.sigma_mean_vcd(coincidences(:,2)).^2 +...
                table2.std_vcd(coincidences(:,2)).^2 );
end

% time differences in hours
time_diff=coincidences(:,3)*24;

% average time of each coincidence
times=(table1.mjd2k(coincidences(:,1))+table2.mjd2k(coincidences(:,2)))/2;


% save sza info
sza1=table1.sza(coincidences(:,1));
sza2=table2.sza(coincidences(:,2));

% get DMP info
T1=NaN(length(sza1),12);
T2=NaN(length(sza2),12);
spv1=NaN(length(sza1),12);
spv2=NaN(length(sza2),12);
lat1=NaN(length(sza1),12);
lat2=NaN(length(sza2),12);

% get profiles
prof1=NaN(length(sza1),38);
prof2=NaN(length(sza2),38);

try T1=table1.T(coincidences(:,1),:); spv1=table1.spv(coincidences(:,1),:); end

try
    lat1=table1.lat(coincidences(:,1),:);
catch
    try lat1=table1.lat_dmp(coincidences(:,1),:); end
end

try prof1=table1.num_dens(coincidences(:,1),:); end
%%%

try T2=table2.T(coincidences(:,2),:); spv2=table2.spv(coincidences(:,2),:); end

try
    lat2=table2.lat(coincidences(:,2),:);
catch
    try lat2=table2.lat_dmp(coincidences(:,2),:); end
end

try prof2=table2.num_dens(coincidences(:,2),:); end


%% Sort results by time
[times,sortind]=sort(times);

vcd1=vcd1(sortind);
vcd2=vcd2(sortind);

error1=error1(sortind);
error2=error2(sortind);

sza1=sza1(sortind);
sza2=sza2(sortind);

T1=T1(sortind,:);
T2=T2(sortind,:);
spv1=spv1(sortind,:);
spv2=spv2(sortind,:);
lat1=lat1(sortind,:);
lat2=lat2(sortind,:);

DMP=struct();
DMP.sza1=sza1;
DMP.sza2=sza2;
DMP.T1=T1;
DMP.T2=T2;
DMP.spv1=spv1;
DMP.spv2=spv2;
DMP.lat1=lat1;
DMP.lat2=lat2;

time_diff=time_diff(sortind);



% debug: show plot of two time series, with lines connecting the selected
% measurement pairs
% % t1=table1.mjd2k(coincidences(:,1));
% % t2=table2.mjd2k(coincidences(:,2));
% % x=ones(size(t1'));
% % 
% % figure
% % plot([t1';t2'],[x;x*-1], 'k.-'), hold on
% % plot(table1.mjd2k', ones(size(table1.mjd2k')),'ro')
% % plot(table2.mjd2k', ones(size(table2.mjd2k'))*-1,'ro')
% % ylim([-5,5])

end