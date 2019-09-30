function [ vcd1, vcd2, error1, error2, times] = find_coincidences_twilight( table1, table2, partcol )
%[ vcd1, vcd2, error1, error2 ] = find_coincidences_twilight( table1, table2, partcol )
%   Find coincident GBS/SAOZ/ACE measurements (same twilight)
%
%   Input needs to be a table in the standard GBS format (first 3 columns
%   are year, doy, ampm). Satellite tables are produced by
%   reformat_ACE_FTS.m, 
%
%   Optional: partcol: if 'true', finds partial column coincidences from
%   satellite occultation measurements

if nargin==2
    partcol=false;
end

% check input tables to see if data is GB or sat
one_doas=false;
two_doas=false;

if ~partcol
    if any(strcmp('mean_vcd',table1.Properties.VariableNames)==1), one_doas=true; end
    if any(strcmp('mean_vcd',table2.Properties.VariableNames)==1), two_doas=true; end
end

%% find coincidences
% coincidences complicated by the fact that ACE might have multiple
% measurements during the same twilight tht are close to PEARL
if one_doas || two_doas

    % find coincidences (first ones in case of multiple meas.)
    arr1=table1(:,1:3);
    arr2=table2(:,1:3);

    [~,ind1,ind2]=intersect(arr1,arr2,'rows');
    
    % find additional satellite measurements
    if ~one_doas
        
        % indices of all repeated meas.
        extra=find_duplicates(table1);
        
        % find whether twilights are actually coincidences
        [~,tmp]=ismember(extra(:,2),ind1);
        
        % remove non-coincidences
        extra(tmp==0,:)=[];
        tmp(tmp==0,:)=[];

        % add table1 indices
        ind1=[ind1; extra(:,1)];
        
        % add table2 indices
        ind2=[ind2; ind2(tmp)];
        
    % same as above, but reversed order
    elseif ~two_doas
        
        extra=find_duplicates(table2);
        
        [~,tmp]=ismember(extra(:,2),ind2);
        
        extra(tmp==0,:)=[];
        tmp(tmp==0,:)=[];

        ind2=[ind2; extra(:,1)];
        
        ind1=[ind1; ind1(tmp)];
        
    end
else

    % find coincidences for FTS and MAESTRO: same time info, so time info
    % is also passed to intersect -- same occultations are matched
    arr1=table1(:,1:4);
    arr2=table2(:,1:4);

    [~,ind1,ind2]=intersect(arr1,arr2,'rows');
    
end



%% allocate variables

% avg time of each coincidence
times=(table1.mjd2k(ind1)+table2.mjd2k(ind2))/2;

% sort indices by time
[times,sortind]=sort(times);

ind1=ind1(sortind);
ind2=ind2(sortind);

% columns and errors
if ~partcol
    if one_doas
        % GBS table
        % VCDs
        vcd1=table1.mean_vcd(ind1);
        % VCD errors
        error1=sqrt(table1.sigma_mean_vcd(ind1).^2 + table1.std_vcd(ind1).^2);

    else    
        % satellite table
        % extended satellite column
        vcd1=table1.tot_col(ind1);
        % error
        error1=table1.tot_col_err(ind1);
    end

    if two_doas
        % GBS table
        % VCDs
        vcd2=table2.mean_vcd(ind2);
        % VCD errors
        error2=sqrt(table2.sigma_mean_vcd(ind2).^2 + table2.std_vcd(ind2).^2);

    else    
        % satellite table
        % extended satellite column
        vcd2=table2.tot_col(ind2);
        % error
        error2=table2.tot_col_err(ind2);
    end
    
else
    
    % extended satellite column
    vcd1=table1.part_col(ind1);
    vcd2=table2.part_col(ind2);
    % error
    error1=table1.part_col_err(ind1);
    error2=table2.part_col_err(ind2);

end
end


%%
function extra_meas = find_duplicates(table_in)
    % find satellite measurements of the same twilight
    
    % loop through all the data
    search=true;
    i=1;
    extra_meas=[0,0];
    
    while search
        
        % skip repetitions that are already found
        if any(extra_meas(:,1)==i), i=i+1; continue, end
        
        % find rows matching current one
        ind_tmp=find(table_in.year==table_in.year(i) & table_in.day==table_in.day(i) &...
                     table_in.ampm==table_in.ampm(i));
        
        if length(ind_tmp)>1 % matches are found
            % find indices of extra measurements excluding current one
            % (intersect will find the first index, just like this code)
            tmp=setdiff(ind_tmp,i);
            
            % save extra indices, and also the index they duplicate
            % (to match with intersect)
            extra_meas=[extra_meas; [tmp, ones(size(tmp))*i]];
        end
        
        % stop if end of table is reached
        if i==size(table_in,1), search=false; end
        
        i=i+1;
        
    end
    
    % remove initial line
    extra_meas(1,:)=[];
    
end
    