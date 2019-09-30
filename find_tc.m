function [ R_t, RMSE, SNR, N ] = find_tc( table1, type1, table2, type2, table3, type3 )
%FIND_TC find triple coincidences based on twilight/time coincidence criteria
%
%   Input: 
%       Tables for three datasets, onse satellite and two ground-based
%       Type of each dataset
%           'tw': twilight DOAS data
%           'ds': direct sun (bruker or brewer) data
%           'ace': ACE-FTS or ACE-MAESTRO, twilight data
%           'os': OSIRIS, 
%       Code selects coincidence function based on type of datasets, using
%       total columns only (!)

%% setup

dt=12;
partcol=false; 

%% find pairwise coincidences based on types of individual datasets

%%% pair 1
cell1={type1, type2};

if find_in_cell(cell1,'tw',true) && (strcmp(type1,type2) || find_in_cell(cell1,'ace',true))
    % two DOAS instruments; or DOAS and ACE
    [ x12, x21, e12, e21 ] = find_coincidences_twilight( table1, table2, partcol );
elseif find_in_cell(cell1,'ace',true) && strcmp(type1,type2)
    % ACE-FTS and ACE-MAESTRO, compare total columns, same occultation only
    [ x12, x21, e12, e21 ] = find_coincidences_twilight( table1, table2, partcol, true );
else
    % time-based coincidences
    [ x12, x21, e12, e21 ] = find_coincidences_time( table1, table2, dt, partcol );
end

%%% pair 2
cell2={type2, type3};

if find_in_cell(cell2,'tw',true) && (strcmp(type2,type3) || find_in_cell(cell2,'ace',true))
    % two DOAS instruments; or DOAS and ACE
    [ x23, x32, e23, e32 ] = find_coincidences_twilight( table2, table3, partcol );
elseif find_in_cell(cell2,'ace',true) && strcmp(type2,type3)
    % ACE-FTS and ACE-MAESTRO, compare total columns, same occultation only
    [ x23, x32, e23, e32 ] = find_coincidences_twilight( table2, table3, partcol, true );
else
    % time-based coincidences
    [ x23, x32, e23, e32 ] = find_coincidences_time( table2, table3, dt, partcol );
end

%%% pair 3
cell3={type1, type3};

if find_in_cell(cell3,'tw',true) && (strcmp(type1,type3) || find_in_cell(cell3,'ace',true))
    % two DOAS instruments; or DOAS and ACE
    [ x13, x31 ] = find_coincidences_twilight( table1, table3, partcol );
elseif find_in_cell(cell3,'ace',true) && strcmp(type1,type3)
    % ACE-FTS and ACE-MAESTRO, compare total columns, same occultation only
    [ x13, x31 ] = find_coincidences_twilight( table1, table3, partcol, true );
else
    % time-based coincidences
    [ x13, x31 ] = find_coincidences_time( table1, table3, dt, partcol );
end

    
%% construct triple coincidence matrix
% use value-error paris to match matrices, just in case 2 values are equal

% start from coincidences between instruments 1 and 2, and check for
% repeated values in coincidences between instruments 2 and 3
[~,ind21,ind23]=intersect_repeat([x21,e21],[x23,e23]);

% array with each column corresponding to one instrument
%   for now, matrix contains measurements when instrument 2 has coincident
%   measurements with both instr1 and instr3
triple=[x12(ind21), x21(ind21), x32(ind23)];
triple_err=[e12(ind21), e21(ind21), e32(ind23)];

% make sure empty matrix doesn't break code (output is NaN)
triple=[[0,0,0]; triple];
triple_err=[[0,0,0]; triple_err]; 

% need to check whether 1-3 pairs in triple are actualy coincidences (will
% be for twilight, but might not be for time)

mask=ismember(triple(:,1:2:3), [x13,x31],'rows');

triple=triple(mask,:);
triple_err=triple_err(mask,:);


%% Calculate statistics

for i=1:3
        
    inds= setdiff([1,2,3],i);

    % list of covariance matrices: current vs other1, current vs other2, other1 vs other2
    tmp1=cov(triple(:,[i,inds(1)]));
    tmp2=cov(triple(:,[i,inds(2)]));
    tmp3=cov(triple(:,inds));
    
    % ratio that appears in TC formulas
    ratio=tmp1(1,2)*tmp2(1,2)/tmp3(1,2);

    % variance of current
    var_i=tmp1(1,1);

    if var_i>ratio
    
        R_t(i)=sqrt(ratio/var_i);

        RMSE(i)=sqrt(var_i - ratio );

        SNR(i)=-10*log10(var_i/ratio - 1 );

    else % something's wrong
        
        eval(['type' num2str(i) '=[type' num2str(i) ' ''*''];'])
        
        warning(['R_t > 1 for starred member in triplet; NaNs returned (' type1 '-' type2 '-' type3 ')']);
        R_t(i)=NaN;
        RMSE(i)=NaN;
        SNR(i)=NaN;
        
        
    end
end

N=size(triple,1);

end

