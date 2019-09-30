function [ diff ] = rel_abs_diff( type, x1, x2 )
%[ diff ] = re_abs_diff( type, x1, x2 ) Calculates absolute or relative difference
%
% INPUT:    type: 'abs' or 'rel', for absolute and relative difference. Relative 
%               difference is taken w.r.t the mean of the two datapoints
%           x1, x2: coincident measurements (must be same size)
%
% OUTPUT:   diff: rel or abs difference between each datapoint
%   
% Kristof Bognar, April 2018

%% consistency checks
if length(x1)~=length(x2), error('Input arrays must be the same size'); end

% if any(isnan(x1+x2)), error('Input data contains NaNs'); end

%% calculate values

switch type
    case 'abs'
        % diff=abs(x1-x2);
        diff=x1-x2;
    case 'rel'
        % diff=abs(x1-x2) ./ ((x1+x2)./2);
        diff=( (x1-x2) ./ ((x1+x2)./2) )*100; % convert to percent
    otherwise
        error('Type option must be either "rel" or "abs"')
end

end

