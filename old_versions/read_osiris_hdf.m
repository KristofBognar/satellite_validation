function [data] = read_osiris_hdf( filename )
%% functions from OSIRIS website (for HDF data)

    swathinfo =  h5info(filename, '/HDFEOS/SWATHS');
    variables = L_LoadGroupVariables( filename, swathinfo);
    
    nvar  = numel(variables);
    data = [];
    if (nvar > 0)
        for i = 1:nvar
            name = variables{i}.name;
            data.(name) = variables{i}.value;
            bad = (data.(name) == -9999.0);
            data.(name)(bad) = NaN;
        end
    end
end
%%%
function [variables] = L_LoadGroupVariables( filename, group)

    variables = {};
    for i = 1: numel(group.Groups)
        [v] = L_LoadGroupVariables( filename, group.Groups(i));   
        variables = [variables; v];
    end
    
    for i = 1: numel(group.Datasets)
        v = L_LoadDataSet( filename, group.Name, group.Datasets(i) );
        variables = [variables; v];
    end
end
%%%
function [v] = L_LoadDataSet( filename, groupname, dataset )

    entry       = [];
    datasetname = [ groupname, '/',dataset.Name];
    value       = h5read( filename, datasetname );
    entry.name  = dataset.Name;
    entry.value = value;
    v           = {entry};
end

