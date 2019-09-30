% function [ vcd1, vcd2, error1, error2, time_diff ] = find_coincidences_time( table1, table2, dt, partcol )
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
%       
%   Note: function eats memory and processing power if data dimensions are
%       in the 10 000s
%
% Kristof Bognar, January 2018

if 0
table1=saoz;
table2=brewer_ds;
else
table1=brewer_ds;
table2=osiris;
end

dt=0.5;
partcol=false;


% max time difference between measurements 
% if nargin==2
%     dt=0.5; % default: 12 hours
%     partcol=false;
% elseif nargin==3
%     dt=dt/24;
%     partcol=false;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% find temporal coincidences

% % % % copy each array into 2D with dimensions of table2 x table1
% % % tmp1=repmat(table1.mjd2k',length(table2.mjd2k),1);
% % % tmp2=repmat(table2.mjd2k,1,length(table1.mjd2k));
% % % % find differences between each element
% % % diff=abs(tmp2-tmp1);

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

% find unique indices (i.e. each element with repetitions removed, NOT
% elements that don't repeat)
% ia is first occurence of each element
[tmp,ia,~]=unique(min_ind12);
% find how many times each index repeats
reps=hist(min_ind12,unique(min_ind12));

% do the same the other way, 

% save indices that only appear once, those are coincidences
% coincidences array has columns:
%   table1 indices
%   table2 indices
%   time difference in hours
coincidences=[ia(reps==1),tmp(reps==1)',min_val12(ia(reps==1))']; 

% deal with repeated indices
% look up what's the closest time to each repeated table2 index -- min_ind_21 gives 
% the index of closest times to each table2 time
coincidences=[coincidences; [min_ind21(tmp(reps>1)),tmp(reps>1)',min_val21(tmp(reps>1))]];

% % discard coincidences outside time window
coincidences(coincidences(:,3)>dt,:)=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %% Assign output values
% % % 
% % % % figure out date type for each table, based on unique columns
% % % % 1: satelite, 2: bruker, 3: brewer 4: GBS
% % % 
% % % if any(strcmp('dist',table1.Properties.VariableNames)==1)
% % %     type1=1;
% % % elseif any(strcmp('tot_col_err_sys',table1.Properties.VariableNames)==1)
% % %     type1=2;
% % % elseif any(strcmp('tot_col_err',table1.Properties.VariableNames)==1)
% % %     type1=3;
% % % elseif any(strcmp('mean_vcd',table1.Properties.VariableNames)==1)
% % %     type1=4;
% % % end
% % % 
% % % if any(strcmp('dist',table2.Properties.VariableNames)==1)
% % %     type2=1;
% % % elseif any(strcmp('tot_col_err_sys',table2.Properties.VariableNames)==1)
% % %     type2=2;
% % % elseif any(strcmp('tot_col_err',table2.Properties.VariableNames)==1)
% % %     type2=3;
% % % elseif any(strcmp('mean_vcd',table2.Properties.VariableNames)==1)
% % %     type2=4;
% % % end
% % % 
% % % % assign VCD and error values
% % % % table 1
% % % if type1==1
% % %     if ~partcol
% % %         vcd1=table1.tot_col(coincidences(:,1));
% % %         error1=table1.tot_col_err(coincidences(:,1));
% % %     else
% % %         vcd1=table1.part_col(coincidences(:,1));
% % %         error1=table1.part_col_err(coincidences(:,1));
% % %     end        
% % % elseif type1==2
% % %     vcd1=table1.tot_col(coincidences(:,1));
% % %     error1=sqrt(table1.tot_col_err_sys(coincidences(:,1)).^2 +...
% % %                 table1.tot_col_err_rand(coincidences(:,1)).^2 );
% % % elseif type1==3
% % %     vcd1=table1.tot_col(coincidences(:,1));
% % %     error1=table1.tot_col_err(coincidences(:,1));
% % % elseif type1==4
% % %     vcd1=table1.mean_vcd(coincidences(:,1));
% % %     error1=sqrt(table1.sigma_mean_vcd(coincidences(:,1)).^2 +...
% % %                 table1.std_vcd(coincidences(:,1)).^2 );
% % % end
% % % 
% % % % table 2
% % % if type2==1
% % %     if ~partcol
% % %         vcd2=table2.tot_col(coincidences(:,2));
% % %         error2=table2.tot_col_err(coincidences(:,2));
% % %     else
% % %         vcd2=table2.part_col(coincidences(:,2));
% % %         error2=table2.part_col_err(coincidences(:,2));
% % %     end        
% % % elseif type2==2
% % %     vcd2=table2.tot_col(coincidences(:,2));
% % %     error2=sqrt(table2.tot_col_err_sys(coincidences(:,2)).^2 +...
% % %                 table2.tot_col_err_rand(coincidences(:,2)).^2 );
% % % elseif type2==3
% % %     vcd2=table2.tot_col(coincidences(:,2));
% % %     error2=table2.tot_col_err(coincidences(:,2));
% % % elseif type2==4
% % %     vcd2=table2.mean_vcd(coincidences(:,2));
% % %     error2=sqrt(table2.sigma_mean_vcd(coincidences(:,2)).^2 +...
% % %                 table2.std_vcd(coincidences(:,2)).^2 );
% % % end
% % % 
% % % % time differences in hours
% % % time_diff=coincidences(:,3)*24;

% end