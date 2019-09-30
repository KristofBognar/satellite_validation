function [rt,rmse,snr,N,mean_all,std_all]=get_tc_stats_o3_no_PARIS(...
          gbs_o3,saoz_o3,bruker_o3,brewer_o3_ds,osiris_o3,ace_fts_o3,ace_mae_o3)
%GET_TC_STATS compile triple collocation statistics for each instrument
%combination

%% setup

if nargin==0
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
    load(data_file)
end

GV=gbs_o3;
SA=saoz_o3;
BK=bruker_o3;
BW=brewer_o3_ds;

OS=osiris_o3;
AF=ace_fts_o3;
AM=ace_mae_o3;

clearvars -except GV SA BK BW OS AF AM

% tables to store data in (keep track of what triplets were used)
% rows: instruments in order described below
% cols: all possible pairs of instruments (data is only saved where relevant)
%       1-6: GB-GB;  7-16: GB-sat;  17-19: sat-sat,
%
% index of each instrument:
%
%               1   2   3   4   5   6   7
%               GV, SA, BK, BW, OS, AF, AM
% twilight:     ^^^^^^              ^^^^^^
% ground-based: ^^^^^^^^^^^^^^            
% satellite:                    ^^^^^^^^^^

header={'GVSA','GVBK','GVBW','SABK','SABW','BKBW',...
        'GVOS','GVAF','GVAM','SAOS','SAAF','SAAM','BKOS','BKAF','BKAM','BWOS',...
        'OSAF','OSAM','AFAM'};
                                       
rt=array2table(NaN(7,19));
rt.Properties.VariableNames=header;                                                            

rmse=rt;
snr=rt;
N=rt;

%% all ground-based
disp('GB triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', BK, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','BK',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', BW, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','BW',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', BW, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','BW',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', BW, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','BW',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

%% 2x GB + 1x sat
disp('2GB+1sat triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( BK, 'ds', BW, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'BK','BW','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% skip BK BW AF
% skip BK BW AM

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BW, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BW','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% skip GV BW AF
% skip GV BW AM

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BW, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BW','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% skip SA BW AF
% skip SA BW AM

%% 1x GB + 2x sat
disp('1GB+2sat triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', OS, 'os', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', OS, 'os', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( BK, 'ds', OS, 'os', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'BK','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', OS, 'os', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','OS','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', OS, 'os', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','OS','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( BK, 'ds', OS, 'os', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'BK','OS','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', AF, 'ace', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','AF','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', AF, 'ace', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','AF','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( BK, 'ds', AF, 'ace', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'BK','AF','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

% skip BW OS AF
% skip BW OS AM
% skip BW AF AM

%% all satellites
disp('Sat triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( OS, 'os', AF, 'ace', AM, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'OS','AF','AM',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

%% average by group type
mean_all=table;

tmp=[nanmean(rt{:,:},2),nanmean(rmse{:,:},2),nanmean(snr{:,:},2)];
mean_all.all=tmp;

tmp=[nanmean(rt{:,1:6},2),nanmean(rmse{:,1:6},2),nanmean(snr{:,1:6},2)];
mean_all.GBx2=tmp;

tmp=[nanmean(rt{:,7:16},2),nanmean(rmse{:,7:16},2),nanmean(snr{:,7:16},2)];
mean_all.GBxSAT=tmp;

tmp=[nanmean(rt{:,17:19},2),nanmean(rmse{:,17:19},2),nanmean(snr{:,17:19},2)];
mean_all.SATx2=tmp;

std_all=table;

tmp=[nanstd(rt{:,:},1,2),nanstd(rmse{:,:},1,2),nanstd(snr{:,:},1,2)];
std_all.all=tmp;

tmp=[nanstd(rt{:,1:6},1,2),nanstd(rmse{:,1:6},1,2),nanstd(snr{:,1:6},1,2)];
std_all.GBx2=tmp;

tmp=[nanstd(rt{:,7:16},1,2),nanstd(rmse{:,7:16},1,2),nanstd(snr{:,7:16},1,2)];
std_all.GBxSAT=tmp;

tmp=[nanstd(rt{:,17:19},1,2),nanstd(rmse{:,17:19},1,2),nanstd(snr{:,17:19},1,2)];
std_all.SATx2=tmp;

end

function [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,instr1,instr2,instr3,rt_tmp,rmse_tmp,snr_tmp,N_tmp)

    instr_all={'GV', 'SA', 'BK', 'BW', 'OS', 'AF', 'AM'};
    
    N_tmp=repmat(N_tmp,1,3);

    for i={'rt','rmse','snr','N'}
        eval([i{1} '.' instr2 instr3 '(' num2str(find_in_cell(instr_all,instr1)) ')=' i{1} '_tmp(1);']); 
        eval([i{1} '.' instr1 instr3 '(' num2str(find_in_cell(instr_all,instr2)) ')=' i{1} '_tmp(2);']); 
        eval([i{1} '.' instr1 instr2 '(' num2str(find_in_cell(instr_all,instr3)) ')=' i{1} '_tmp(3);']); 
    end
    
end

