function reorganize_files( instr ) 

% save current directory
cur_dir=pwd;

switch instr
    case 'osiris'
        data_dir='/home/kristof/work/satellite_validation/ODIN-OSIRIS_data';


        % change to OSIRIS data directory
        cd(data_dir)

        % check if 'daily' folder exists
        if exist('daily','dir')
            work_dir=[data_dir '/daily'];
        else
            work_dir=data_dir;
        end

        cd(work_dir)

        % make list of monthly directories
        dirlist=dir_list();

        % make new folders
        cd(data_dir)
        if ~exist('NO2','dir'), mkdir('NO2'), end
        if ~exist('O3','dir'), mkdir('O3'), end
        cd(work_dir)

        % loop through monthly folders
        for i=1:length(dirlist)
            cd([work_dir '/' dirlist{i}]);

            try movefile('*O3*',[data_dir '/O3']); end
            try movefile('*NO2*',[data_dir '/NO2']); end

        end

    case 'maestro'

        data_dir='/home/kristof/work/satellite_validation/ACE-MAESTRO_data/';
        work_dir=[data_dir 'data_tmp/'];

        % save current directory
        cur_dir=pwd;

        % make new folders
        cd(data_dir)
        if ~exist('NO2','dir'), mkdir('NO2'), end
        if ~exist('O3','dir'), mkdir('O3'), end
        if ~exist('O3UV','dir'), mkdir('O3UV'), end

        % change to directory with the data
        cd(work_dir)

        % make list of monthly directories
        yearlist=dir_list();
        
        % loop through yearly folders
        for yy=1:length(yearlist)

            cd([work_dir yearlist{yy}])
            
            % make list of monthly folders
            monthlist=dir_list();
            
            % loop through monthly folders
            for i=1:length(monthlist)
                cd(monthlist{i});

                try movefile('*vo3g*',[data_dir 'O3']); end
                try movefile('*uo3g*',[data_dir 'O3UV']); end
                try movefile('*uno2g*',[data_dir 'NO2']); end

                cd('../')
            end
        end
end

cd(cur_dir)
end
%%
function [names]=dir_list()
        tmp=dir('*');
        names={tmp.name};
        names(1:2)=[];
end