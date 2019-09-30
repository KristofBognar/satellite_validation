function [ vcd1, vcd2, error1, error2, times, DMP, prof1, prof2] =...
         find_coincidences_twilight( table1, table2, partcol, ace_totcol )
%[ vcd1, vcd2, error1, error2 ] = find_coincidences_twilight( table1, table2, partcol )
%   Find coincident GBS/SAOZ/ACE measurements (same twilight)
%
%   Input needs to be a table in the standard GBS format (first 3 columns
%   are year, doy, ampm). Satellite tables are produced by
%   reformat_ACE_FTS.m, 
%
%   Optional:   partcol: if true, finds partial column coincidences from
%   ACE occultation measurements (same occultation only)
%               ace_totcol: if true, finds total column coincidences from
%   ACE occultation measurements (same occultation only)
%
% Kristof Bognar, December 2017

if nargin<4
    ace_totcol=false;
    if nargin==2
        partcol=false;
    end
end

if ace_totcol, partcol=false; end

% check input tables to see if data is GB or sat
one_doas=false;
two_doas=false;

if ~partcol
    if any(strcmp('mean_vcd',table1.Properties.VariableNames)==1), one_doas=true; end
    if any(strcmp('mean_vcd',table2.Properties.VariableNames)==1), two_doas=true; end
end

%% find coincidences
if partcol || ace_totcol % ACE-FTS vs ACE-MAESTRO, compare same occultation only

    % set up arrays (use year, day, ampm and fractional day fields --
    % select specific measurements)
    arr1=table1(:,1:4);
    arr2=table2(:,1:4);

    % find coincidences (each line that appears in both arrays -- each
    % array contains unique lines only)
    [~,ind1,ind2]=intersect(arr1,arr2);
    
else % ACE vs ground-based DOAS
    % coincidences complicated by the fact that ACE might have multiple
    % measurements during the same twilight that are close to PEARL -- use my
    % intersect function that keeps repetitions

    % set up arrays (use year, day and ampm fields only -- select twilights)
    arr1=table1(:,1:3);
    arr2=table2(:,1:3);

    % find coincidences (non-unique twilights are paired in all
    % possible permutations)
    [~,ind1,ind2]=intersect_repeat(arr1,arr2);

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
    
    % partial satellite column
    vcd1=table1.part_col(ind1);
    vcd2=table2.part_col(ind2);
    % error
    error1=table1.part_col_err(ind1);
    error2=table2.part_col_err(ind2);

end

% save sza info
sza1=table1.sza(ind1);
sza2=table2.sza(ind2);

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


try T1=table1.T(ind1,:); spv1=table1.spv(ind1,:); end

% do tangent lat for satellites, dmp derived lat otherwise
try
    lat1=table1.lat(ind1,:);
catch
    try lat1=table1.lat_dmp(ind1,:); end
end

try prof1=table1.num_dens(ind1,:); end

try T2=table2.T(ind2,:); spv2=table2.spv(ind2,:); end

try
    lat2=table2.lat(ind2,:);
catch
    try lat2=table2.lat_dmp(ind2,:); end
end

try prof2=table2.num_dens(ind2,:); end

DMP=struct();
DMP.sza1=sza1;
DMP.sza2=sza2;
DMP.T1=T1;
DMP.T2=T2;
DMP.spv1=spv1;
DMP.spv2=spv2;
DMP.lat1=lat1;
DMP.lat2=lat2;

end


