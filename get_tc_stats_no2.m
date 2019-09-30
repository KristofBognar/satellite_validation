function [rt,rmse,snr,N,mean_all,std_all]=get_tc_stats_no2(...
          gbs_no2,gbs_no2uv,saoz_o3,bruker_no2,osiris_no2,ace_fts_no2)
%GET_TC_STATS compile triple collocation statistics for each instrument
%combination

%% setup

if nargin==0
    data_file='/home/kristof/work/satellite_validation/all_data_nosmooth.mat';
    load(data_file)
end

GV=gbs_no2;
GU=gbs_no2uv;
SA=saoz_o3;
BK=bruker_no2;
BK.tot_col=BK.part_col; % use 12-40 km column for comparisons

OS=osiris_no2;
AF=ace_fts_no2;

clearvars -except GV GU SA BK OS AF

% tables to store data in (keep track of what triplets were used)
% rows: instruments in order described below
% cols: all possible pairs of instruments (data is only saved where relevant)
%       1-6: GB-GB;  7-16: GB-sat;  17-19: sat-sat
%
% index of each instrument:
%
%               1   2   4   3   5   6 
%               GV, GU, SA, BK, OS, AF
% twilight:     ^^^^^^^^^^          ^^
% ground-based: ^^^^^^^^^^^^^^
% satellite:                    ^^^^^^

header={'GVGU','GVSA','GVBK','GUSA','GUBK','SABK',...
        'GVOS','GVAF','GUOS','GUAF','SAOS','SAAF','BKOS','BKAF',...
        'OSAF'};
                                       
rt=array2table(NaN(6,15));
rt.Properties.VariableNames=header;                                                            

rmse=rt;
snr=rt;
N=rt;

%% all ground-based
disp('GB triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', GU, 'tw', SA, 'tw' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','GU','SA',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', GU, 'tw', BK, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','GU','BK',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', BK, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','BK',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
    
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', SA, 'tw', BK, 'ds' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','SA','BK',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

    
%% 2x GB + 1x sat
disp('2GB+1sat triplets')

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', GU, 'tw', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','GU','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', GU, 'tw', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','GU','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', SA, 'tw', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','SA','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', BK, 'ds', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','BK','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', SA, 'tw', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','SA','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', SA, 'tw', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','SA','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', BK, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','BK','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', BK, 'ds', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','BK','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
    
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', OS, 'os' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','OS',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
[rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', BK, 'ds', AF, 'ace' );
    [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','BK','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

% % not enough OS-AF coincidences for this    
% % % %% 1x GB + 2x sat
% % % disp('1GB+2sat triplets')
% % % 
% % % [rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GV, 'tw', OS, 'os', AF, 'ace' );
% % %     [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GV','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% % % [rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( GU, 'tw', OS, 'os', AF, 'ace' );
% % %     [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'GU','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% % % [rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( SA, 'tw', OS, 'os', AF, 'ace' );
% % %     [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'SA','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);
% % % % [rt_tmp,rmse_tmp,snr_tmp,N_tmp]=find_tc( BK, 'ds', OS, 'os', AF, 'ace' );
% % % %     [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,'BK','OS','AF',rt_tmp,rmse_tmp,snr_tmp,N_tmp);

    
    
%% average by group type
mean_all=table;

tmp=[nanmean(rt{:,:},2),nanmean(rmse{:,:},2),nanmean(snr{:,:},2)];
mean_all.all=tmp;

tmp=[nanmean(rt{:,1:6},2),nanmean(rmse{:,1:6},2),nanmean(snr{:,1:6},2)];
mean_all.GBx2=tmp;

tmp=[nanmean(rt{:,7:14},2),nanmean(rmse{:,7:14},2),nanmean(snr{:,7:14},2)];
mean_all.GBxSAT=tmp;

tmp=[nanmean(rt{:,15},2),nanmean(rmse{:,15},2),nanmean(snr{:,15},2)];
mean_all.SATx2=tmp;

std_all=table;

tmp=[nanstd(rt{:,:},1,2),nanstd(rmse{:,:},1,2),nanstd(snr{:,:},1,2)];
std_all.all=tmp;

tmp=[nanstd(rt{:,1:6},1,2),nanstd(rmse{:,1:6},1,2),nanstd(snr{:,1:6},1,2)];
std_all.GBx2=tmp;

tmp=[nanstd(rt{:,7:14},1,2),nanstd(rmse{:,7:14},1,2),nanstd(snr{:,7:14},1,2)];
std_all.GBxSAT=tmp;

tmp=[nanstd(rt{:,15},1,2),nanstd(rmse{:,15},1,2),nanstd(snr{:,15},1,2)];
std_all.SATx2=tmp;

end

function [rt,rmse,snr,N]=fill_tables(rt,rmse,snr,N,instr1,instr2,instr3,rt_tmp,rmse_tmp,snr_tmp,N_tmp)

    instr_all={'GV', 'GU', 'SA', 'BK', 'OS', 'AF'};
    
    N_tmp=repmat(N_tmp,1,3);

    for i={'rt','rmse','snr','N'}
        eval([i{1} '.' instr2 instr3 '(' num2str(find_in_cell(instr_all,instr1)) ')=' i{1} '_tmp(1);']); 
        eval([i{1} '.' instr1 instr3 '(' num2str(find_in_cell(instr_all,instr2)) ')=' i{1} '_tmp(2);']); 
        eval([i{1} '.' instr1 instr2 '(' num2str(find_in_cell(instr_all,instr3)) ')=' i{1} '_tmp(3);']); 
    end
    
end

